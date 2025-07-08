'''Goals (class Corrector):
- implements corrections for distances matrices
- class variable for original distanceMatrix
- class variable to store each corrected matrix (one for each correction)
- functions for each correction
'''
import math
import numpy as np 
from scipy.optimize import minimize
import core.PerameterOptimizer as PerameterOptimizer

class Corrector:

    def __init__(self, distanceMatrix, sequenceFile = None):
        self.originalDistances = distanceMatrix
        self.JukesCantorMatrix = []
        self.F81Matrix = []
        self.K80Matrix = []
        self.T92Matrix = []
        self.TN93Matrix = []
        self.GTRMatrix = []


        #Needed for F81 and F81 MOD
        #Counts frequencies for each of the nucleotides if the file exists
        if sequenceFile is not None:
            self.F81Matrix = []
            self.F81ModMatrix = []
            seqFile = open(sequenceFile)
            content2 = seqFile.readlines()
            self.sequences = content2
            self.frequenciesForFile = {"A" : 0, "T": 0, "G": 0, "C": 0, "Total": 0}
            self.frequenciesPerSequence = [] #list of dicts: each inner dict i contains the ATGC freqs for the associated sequence i
            self.findFrequencies()
            
            


    def findFrequencies(self):
        # We need total nucleotides and totalA,T,G,C for the entire sequence file for F81
        # We need total nucleos and total A,T,G,C for each individual sequence as well for F81 Mod
        # ASK: what is N in the unaligned sequence? do we include that in the seq Length??

        for s in range(1, len(self.sequences), 2):
            #seqLength = len(self.sequences[s])
            #keeps track of totals for each individual sequence (could be put in dictionary to be cleaner)
            freqsForSequence = {"A" : 0, "T": 0, "G": 0, "C": 0, "Total": 0}
            # totalNucleo = 0
            # freqA = 0
            # freqT = 0
            # freqG = 0
            # freqC = 0
            for ch in self.sequences[s]:
                if ch in ("A","T","G","C"):
                    self.frequenciesForFile[ch] += 1
                    self.frequenciesForFile["Total"] += 1
                    freqsForSequence[ch] += 1
                    freqsForSequence["Total"] += 1
            self.frequenciesPerSequence.append(freqsForSequence)
            
    

    def jukesCantor(self):
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        jukesOriginal = np.array(data, dtype=float)

        log_func = 1 - ((4/3) * jukesOriginal)
        log_func = np.where(log_func <= 0, 0.01, log_func)  # Avoid log(0) or log(negative)

        corrected = -0.75 * np.log(log_func)
        self.JukesCantorMatrix = [[label] + list(row) for label, row in zip(labels, corrected)]
        return self.JukesCantorMatrix
        
        
    def K80(self):
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        K80Original = np.array(data, dtype=float)
        
        guess_p = 0.15
        guess_q = 0.01
        
        print("Initial P: " + str(guess_p) + " Initial Q: " + str(guess_q))
        
        initialPremeters = [guess_p, guess_q]
        optimizer = PerameterOptimizer.PeramterOptimizer(originalDistances=self.originalDistances, sequenceFile=self.frequenciesForFile)
        bounds = [(0, None), (0, None)]
        optimizedPerameters = optimizer.optimize(initialPremeters, modelNumber=2, perameterBounds=bounds)
        
        
        optimized_p = optimizedPerameters[0]
        optimized_q = optimizedPerameters[1]
        
        print("Optimized P: " + str(optimized_p) + " Optimized Q: " + str(optimized_q))
        
        firstLog = abs(1 - 2*optimized_p - optimized_q)
        root = abs(1 - 2*optimized_q)
        if firstLog == 0:
            firstLog = 0.01
        if root == 0:
            root = 0.01
        
        #Apply the K80 correction with optimized perameters
        correctedK80Matrix = K80Original * (-0.50 * math.log(firstLog) * math.sqrt(root))
        
        self.K80Matrix = [[label] + list(row) for label, row in zip(labels, correctedK80Matrix)]
        
        return self.K80Matrix

    
    def K81(self):
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        K81Original = np.array(data, dtype=float)
        
        alpha = 0.15
        beta = 0.01
        gamma = 0.01
        
        print("Initial Alpha: " + str(alpha) + " Initial Beta: " + str(beta) + " Initial Gamma: " + str(gamma))
        
        initialPremeters = [alpha, beta, gamma]
        optimizer = PerameterOptimizer.PeramterOptimizer(originalDistances=self.originalDistances, sequenceFile=self.frequenciesForFile)
        bounds = [(0, None), (0, None), (0, None)]
        
        optimizedPerameters = optimizer.optimize(modelNumber=3, initialPerams=initialPremeters, perameterBounds=bounds)
        
        alpha = optimizedPerameters[0]
        beta = optimizedPerameters[1]
        gamma = optimizedPerameters[2]
        
        print("Optimized Alpha: " + str(alpha) + " Optimized Beta: " + str(beta) + " Optimized Gamma: " + str(gamma))
        
        #Correct the K81 matrix with optimized perameters
        lamb_1 = -4 * (alpha + beta)
        lamb_2 = -4 * (alpha + gamma)
        lamb_3 = -4 * (beta + gamma)
        
        P = 1 - math.exp(lamb_1) - math.exp(lamb_2) - math.exp(lamb_3)
        Q = 1 - math.exp(lamb_1) + math.exp(lamb_2) - math.exp(lamb_3)
        R = 1 - math.exp(lamb_1) + math.exp(lamb_2) - math.exp(lamb_3)
        
        log_func = abs((1 - 2*P - 2*Q) * (1 - 2*P - 2*R) * (1 - 2*Q - 2*R))
        
        if log_func == 0:
            log_func = 0.01
            print("Log func = 0")
        
        correctedK81Matrix = K81Original * (-(1/4) * math.log(log_func))
        
        self.K81Matrix = [[label] + list(row) for label, row in zip(labels, correctedK81Matrix)]
        
        return self.K81Matrix
    
    def T92(self):
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        T92InitialMatrix = np.array(data, dtype=float)
        
        
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        G_and_C = (freqG + freqC)/2
        #Where p is the proportion of sites that show transitional differences and q is the proportion of sites that show transversional differences.
        
        #NOTE: Transitions are 15x more likly than transversions
        p_guess = 0.15
        q_guess = 0.01
        
        print("Initial P: " + str(p_guess) + " Initial Q: " + str(q_guess))
        
        #3 total perameters needed
        initialPremeters = [p_guess, q_guess, G_and_C]
        
        #Returns list of optimized perameters
        optimizer = PerameterOptimizer.PeramterOptimizer(originalDistances=self.originalDistances, sequenceFile=self.frequenciesForFile)
        bounds = [(0, None), (0, None), (0, None)] 
        optimizedPerameters = optimizer.optimize(modelNumber=4, initialPerams=initialPremeters, perameterBounds=bounds)
        
        #Optimized perameters to be used in final matrix
        p_optimized = optimizedPerameters[0]
        q_optimized = optimizedPerameters[1]
        
        print("Optimized P: " + str(p_optimized) + " Optimized Q: " + str(q_optimized))
        
        #Ensuring the log values are positive values
        h = 2 * G_and_C * (1 - G_and_C)
        firstLog = abs(1 - (p_optimized/h) - q_optimized)
        secondLog = abs(1 - (2*q_optimized))
        if firstLog == 0:
            firstLog = 0.01
        if secondLog == 0:
            secondLog = 0.01
            
        d = -h * math.log(firstLog) - 1/2 * (1 - h) * math.log(secondLog)
        
        print("First Log:" + str(firstLog))
        print("Second Log:" + str(secondLog))
        
        correctedT92Matrix = T92InitialMatrix * d
        self.T92Matrix = [[label] + list(row) for label, row in zip(labels, correctedT92Matrix)]
        
        return self.T92Matrix
        
    
    def TN93(self):
        labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        TN93InitialMatrix = np.array(data, dtype=float)
        
        freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
        freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        
        print("Frequencies: A: " + str(freqA) + " T: " + str(freqT) + " G: " + str(freqG) + " C: " + str(freqC))
        
        #Transitions between purines
        alpha_1_initial = 0.01
        #Transitions between pyrimidines
        alpha_2_intial = 0.01
        #Transversions
        beta_initial = 0.15
        #Mean and variance of the number of nucleotide substitutions per site
        a_initial = 0.01
        
        print("Initial Alpha 1: " + str(alpha_1_initial) + " Initial Alpha 2: " + str(alpha_2_intial) + " Initial Beta: " + str(beta_initial) + " Initial A: " + str(a_initial))
        
        initialPerameters = [alpha_1_initial, alpha_2_intial, beta_initial, a_initial]
        
        optimized = PerameterOptimizer.PeramterOptimizer(originalDistances=self.originalDistances, sequenceFile=self.frequenciesForFile)
        bounds = [(0, None), (0, None), (0, None), (0, None)]
        optimizedPerameters = optimized.optimize(modelNumber=5, initialPerams=initialPerameters, perameterBounds=bounds)
        
        alpha_1_optimized = optimizedPerameters[0]
        alpha_2_optimized = optimizedPerameters[1]
        beta_optimized = optimizedPerameters[2]
        a_optimized = optimizedPerameters[3]
        
        print("Optimized Alpha 1: " + str(alpha_1_optimized) + " Optimized Alpha 2: " + str(alpha_2_optimized) + " Optimized Beta: " + str(beta_optimized) + " Optimized A: " + str(a_optimized))
        #Using optimized perameters to set up formulas
        #NOTE: t = 1 for evolutionary time thus not needed here
        
        epsilon = 1e-8  # Small value to avoid division by zero or math errors

        g_r = freqG + freqA
        if abs(g_r) < epsilon:
            print("Warning: g_r is zero or near zero, setting to epsilon.")
            g_r = epsilon

        g_y = freqT + freqC
        if abs(g_y) < epsilon:
            print("Warning: g_y is zero or near zero, setting to epsilon.")
            g_y = epsilon

        Q_estimate = abs(2 * g_r * g_y * (1 - (a_optimized/(a_optimized+2*beta_optimized))**a_optimized))

        # freqA * freqG denominator
        fg_fg = freqA * freqG
        if abs(fg_fg) < epsilon:
            print("Warning: freqA * freqG is zero or near zero, setting to epsilon.")
            fg_fg = epsilon

        # freqT * freqC denominator
        ft_fc = freqT * freqC
        if abs(ft_fc) < epsilon:
            print("Warning: freqT * freqC is zero or near zero, setting to epsilon.")
            ft_fc = epsilon

        # P_1 and P_2 denominators
        #P_1 = abs((2 * freqA * freqG) / g_r * (g_r + g_y * math.exp(-2 * beta_optimized) - math.exp(-2 * (g_r * alpha_1_optimized + g_y * beta_optimized))))
        #P_2 = abs((2 * freqT * freqC) / g_y * (g_y + g_r * math.exp(-2 * beta_optimized) - math.exp(-2 * (g_y * alpha_2_optimized + g_r * beta_optimized))))
        #P_1 = abs(((2 * freqA * freqG) / g_r) * (g_r - (a_optimized / (a_optimized + 2 * (g_r * alpha_1_optimized + g_y * beta_optimized)))**a_optimized) + g_y * (a_optimized/ a_optimized + 2 * beta_optimized)**a_optimized)
        
        #P_2 = abs(((2 * freqT * freqC) / g_y) * (g_y - (a_optimized / (a_optimized + 2 * (g_y * alpha_2_optimized + g_r * beta_optimized)))**a_optimized) + g_r * (a_optimized/ a_optimized + 2 * beta_optimized)**a_optimized)
        epsilon = 1e-8

        # Safe denominator
        g_r_safe = g_r if abs(g_r) > epsilon else epsilon
        g_y_safe = g_y if abs(g_y) > epsilon else epsilon

        # Safe exponentials
        try:
            exp1 = math.exp(-2 * beta_optimized)
        except OverflowError:
            exp1 = 0.0
        try:
            exp2 = math.exp(-2 * (g_r * alpha_1_optimized + g_y * beta_optimized))
        except OverflowError:
            exp2 = 0.0

        P_1_base = g_r + g_y * exp1 - exp2
        if not np.isfinite(P_1_base):
            print(f"Warning: P_1_base is not finite ({P_1_base}), setting to epsilon.")
            P_1_base = epsilon

        P_1 = (2 * freqA * freqG) / g_r_safe * P_1_base

        # Repeat for P_2
        try:
            exp3 = math.exp(-2 * (g_y * alpha_2_optimized + g_r * beta_optimized))
        except OverflowError:
            exp3 = 0.0

        P_2_base = g_y + g_r * exp1 - exp3
        if not np.isfinite(P_2_base):
            print(f"Warning: P_2_base is not finite ({P_2_base}), setting to epsilon.")
            P_2_base = epsilon

        P_2 = (2 * freqT * freqC) / g_y_safe * P_2_base
        #print("P_1: " + str(P_1) + " P_2: " + str(P_2))
        
        # equation_1
        denom1 = abs(2 * fg_fg)
        if abs(denom1) < epsilon:
            print("Warning: 2 * freqA * freqG is zero or near zero, setting to epsilon.")
            denom1 = epsilon
        base1 = abs(1 - (g_r / denom1 * P_1) - (1 / (2 * g_r) * Q_estimate))
        if base1 <= 0:
            print(f"Warning: base1 for equation_1 is <= 0 ({base1}), setting to epsilon.")
            base1 = epsilon
        equation_1 = (freqA * freqG) / g_r * (base1) ** (-1 / a_optimized)

        # equation_2
        denom2 = abs(2 * ft_fc)
        if abs(denom2) < epsilon:
            print("Warning: 2 * freqT * freqC is zero or near zero, setting to epsilon.")
            denom2 = epsilon
        base2 = abs(1 - (g_y / denom2 * P_2) - (1 / (2 * g_y) * Q_estimate))
        if base2 <= 0:
            print(f"Warning: base2 for equation_2 is <= 0 ({base2}), setting to epsilon.")
            base2 = epsilon
        equation_2 = (freqT * freqC) / g_y * (base2) ** (-1 / a_optimized)

        # equation_3
        denom3 = abs(2 * g_r * g_y)
        if abs(denom3) < epsilon:
            print("Warning: 2 * g_r * g_y is zero or near zero, setting to epsilon.")
            denom3 = epsilon
        base3 = abs(1 - (1 / denom3) * Q_estimate)
        if base3 <= 0:
            print(f"Warning: base3 for equation_3 is <= 0 ({base3}), setting to epsilon.")
            base3 = epsilon

        # equation_3 numerator/denominator checks
        gr_denom = abs(g_r)
        if abs(gr_denom) < epsilon:
            print("Warning: g_r in equation_3 denominator is zero or near zero, setting to epsilon.")
            gr_denom = epsilon
        gy_denom = abs(g_y)
        if abs(gy_denom) < epsilon:
            print("Warning: g_y in equation_3 denominator is zero or near zero, setting to epsilon.")
            gy_denom = epsilon

        equation_3 = (g_r * g_y) - ((freqA * freqG * g_y) / gr_denom - (freqT * freqC * g_r) / gy_denom) * (base3) ** (-1 / a_optimized)

        correction = TN93InitialMatrix * abs((2 * a_optimized * (equation_1 + equation_2 + equation_3 - freqA * freqG - freqT * freqC - g_r * g_y)))

        self.TN93Matrix = [[label] + list(row) for label, row in zip(labels, correction)]
        
        return self.TN93Matrix