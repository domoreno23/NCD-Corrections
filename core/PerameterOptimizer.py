'''Interface for optimizing perameters of an evolutionary model
1. Uses arbitrary guesses below 1 for intial perameter values
2. Utilizes a Maximum Likelihood calculation to optimize perameters
3. ONLY returns optimized perameters to be used by the Corrector class
'''
import numpy as np
import math
from scipy.optimize import minimize

class PeramterOptimizer:
    def __init__(self, originalDistances, sequenceFile,):
        #These two will be passed in from the Corrector class
        self.originalDistances = originalDistances
        self.frequenciesForFile = sequenceFile
        self.optimizedPerameters = []
        self.K80Matrix = []
        self.T92Matrix = []
    
    def K80(self, perams):
        data = [row[1:] for row in self.originalDistances]
        K80InitialMatrix = np.array(data, dtype=float)
        p_initial, q_initial = perams[0], perams[1]
        
        firstLog = abs(1 - 2*p_initial - q_initial)
        root = abs(1 - 2*q_initial)
        if firstLog == 0:
            firstLog = 0.01
        if root == 0:
            root = 0.01
        
        self.K80Matrix = K80InitialMatrix * (-0.50 * math.log((firstLog) * math.sqrt(root)))
        
        return np.array(self.K80Matrix)
    
    def K81(self, perams):
        data = [row[1:] for row in self.originalDistances]
        K81InitialMatrix = np.array(data, dtype=float)
        alpha, beta, gamma = perams[0], perams[1], perams[2]
        
        lamb_1 = -4 * (alpha + beta)
        lamb_2 = -4 * (alpha + gamma)
        lamb_3 = -4 * (beta + gamma)
        
        P = (1 - math.exp(lamb_1) - math.exp(lamb_2) + math.exp(lamb_3)) / 4 
        Q = (1 - math.exp(lamb_1) + math.exp(lamb_2) - math.exp(lamb_3)) / 4
        R = (1 + math.exp(lamb_1) - math.exp(lamb_2) - math.exp(lamb_3)) / 4
        
        log_func = abs((1 - 2*P - 2*Q) * (1 - 2*P - 2*R) * (1 - 2*Q - 2*R))
        
        if log_func == 0:
            log_func = 0.01
            print("Log func = 0")
        
        self.K81Matrix = K81InitialMatrix * (-(1/4) * math.log(log_func))
        
        
        return np.array(self.K81Matrix)
    
        
    def T92(self, perams):
        data = [row[1:] for row in self.originalDistances]
        T92InitialMatrix = np.array(data, dtype=float)
        
        p_initial, q_initial, G_and_C = perams[0], perams[1], perams[2]
        
        #Ensuring the log values are not negative
        h = 2 * G_and_C * (1 - G_and_C)
        firstLog = abs(1 - (p_initial/h) - q_initial)
        secondLog = abs(1 - (2*q_initial))
        if firstLog == 0:
            firstLog = 0.01
        if secondLog == 0:
            secondLog = 0.01
            
        d = -h * math.log(firstLog) - 1/2 * (1 - h) * math.log(secondLog)
        
        self.T92Matrix = T92InitialMatrix * d
        
        return np.array(self.T92Matrix)
    
    def TN93(self, perams):
        
        data = [row[1:] for row in self.originalDistances]
        TN93InitialMatrix = np.array(data, dtype=float)
        
        freqA = self.frequenciesForFile["A"] / self.frequenciesForFile["Total"]
        freqT = self.frequenciesForFile["T"] / self.frequenciesForFile["Total"]
        freqG = self.frequenciesForFile["G"] / self.frequenciesForFile["Total"]
        freqC = self.frequenciesForFile["C"] / self.frequenciesForFile["Total"]
        
        alpha_1, alpha_2, beta, a = perams[0], perams[1], perams[2], perams[3]
        epsilon = 1e-8  # Small value to avoid division by zero or math errors

        g_r = freqG + freqA
        if abs(g_r) < epsilon:
            print("Warning: g_r is zero or near zero, setting to epsilon.")
            g_r = epsilon

        g_y = freqT + freqC
        if abs(g_y) < epsilon:
            print("Warning: g_y is zero or near zero, setting to epsilon.")
            g_y = epsilon

        Q_estimate = abs(2 * g_r * g_y * (1 - (a/(a+2*beta))**a))

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
        #P_1 = abs((2 * freqA * freqG) / g_r * (g_r + g_y * math.exp(-2 * beta) - math.exp(-2 * (g_r * alpha_1 + g_y * beta))))
        #P_2 = abs((2 * freqT * freqC) / g_y * (g_y + g_r * math.exp(-2 * beta) - math.exp(-2 * (g_y * alpha_2 + g_r * beta))))
        
        epsilon = 1e-8

        # Safe denominator
        g_r_safe = g_r if abs(g_r) > epsilon else epsilon
        g_y_safe = g_y if abs(g_y) > epsilon else epsilon

        # Safe exponentials
        try:
            exp1 = math.exp(-2 * beta)
        except OverflowError:
            exp1 = 0.0
        try:
            exp2 = math.exp(-2 * (g_r * alpha_1 + g_y * beta))
        except OverflowError:
            exp2 = 0.0

        P_1_base = g_r + g_y * exp1 - exp2
        if not np.isfinite(P_1_base):
            print(f"Warning: P_1_base is not finite ({P_1_base}), setting to epsilon.")
            P_1_base = epsilon

        P_1 = (2 * freqA * freqG) / g_r_safe * P_1_base

        # Repeat for P_2
        try:
            exp3 = math.exp(-2 * (g_y * alpha_2 + g_r * beta))
        except OverflowError:
            exp3 = 0.0

        P_2_base = g_y + g_r * exp1 - exp3
        if not np.isfinite(P_2_base):
            print(f"Warning: P_2_base is not finite ({P_2_base}), setting to epsilon.")
            P_2_base = epsilon

        P_2 = (2 * freqT * freqC) / g_y_safe * P_2_base
        
        #print(f"P_1: {P_1}, P_2: {P_2}")

        # equation_1
        denom1 = abs(2 * fg_fg)
        if abs(denom1) < epsilon:
            print("Warning: 2 * freqA * freqG is zero or near zero, setting to epsilon.")
            denom1 = epsilon
        base1 = abs(1 - (g_r / denom1 * P_1) - (1 / (2 * g_r) * Q_estimate))
        if base1 <= 0:
            print(f"Warning: base1 for equation_1 is <= 0 ({base1}), setting to epsilon.")
            base1 = epsilon
        equation_1 = (freqA * freqG) / g_r * (base1) ** (-1 / a)

        # equation_2
        denom2 = abs(2 * ft_fc)
        if abs(denom2) < epsilon:
            print("Warning: 2 * freqT * freqC is zero or near zero, setting to epsilon.")
            denom2 = epsilon
        base2 = abs(1 - (g_y / denom2 * P_2) - (1 / (2 * g_y) * Q_estimate))
        if base2 <= 0:
            print(f"Warning: base2 for equation_2 is <= 0 ({base2}), setting to epsilon.")
            base2 = epsilon
        equation_2 = (freqT * freqC) / g_y * (base2) ** (-1 / a)

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

        equation_3 = (g_r * g_y) - ((freqA * freqG * g_y) / gr_denom - (freqT * freqC * g_r) / gy_denom) * (base3) ** (-1 / a)

        self.TN93Matrix = TN93InitialMatrix * abs((2 * a * (equation_1 + equation_2 + equation_3 - freqA * freqG - freqT * freqC - g_r * g_y)))
        return np.array(self.TN93Matrix)
        
    
    
    #Calculates differences between model and observed distances, core of MLE
    def loss(self, modelNumber, perams):
        #NOTE: perams are passed through here first and then passed again to each of the model functions
        model_distances = np.array([])
        #labels = [row[0] for row in self.originalDistances]
        data = [row[1:] for row in self.originalDistances]
        observed_distances = np.array(data, dtype=float)
        if modelNumber == 2:
            model_distances = self.K80(perams)
        elif modelNumber == 3:
            model_distances = self.K81(perams)
        elif modelNumber == 4:
            model_distances = self.T92(perams)
        elif modelNumber == 5:
            model_distances = self.TN93(perams)
        
            
        
        
        return np.sum((observed_distances - model_distances) ** 2)
    
    
    #Optimizes perameters by utilizing the differences and initial values
    def optimize(self, initialPerams, modelNumber, perameterBounds):
        initialValues = np.array(initialPerams)
        result = minimize(lambda perams: self.loss(modelNumber, perams), x0=initialValues, bounds=perameterBounds)
        
        optimizedPerams = result.x
        return optimizedPerams