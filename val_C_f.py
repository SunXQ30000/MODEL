# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize

# Parameters
m = 0.6121

Fp = 0.49
Ff = 0.09
Fc = 0.08

b = Fp + Ff + Fc
  
a = 0.2441

Topt = 30
Tmin = 15.8
Tmax = 46.8
Tact = 26

dt = 1  # time step in days

W_observed = np.array([80,130.89736,174.4584,204.76174,334.60311,431])  
t_observed = np.array([0,21,43,58,84,115])  

W0 = W_observed[0]  # 
days = t_observed[-1]  #
W = np.zeros(days + 1)  #
W[0] = W0
dw_dt = np.zeros(days + 1) 
feeding_rates = np.zeros(days + 1)  

def objective(f_values):
    W_temp = np.zeros_like(W)
    W_temp[0] = W0
    for t in range(1, len(W_temp)):
        f = f_values[t-1] 
        dw_dt = f * a * b * W_temp[t-1]**m
        W_temp[t] = W_temp[t-1] + dw_dt * dt
    return np.sum((W_temp[t_observed] - W_observed)**2)  

# 
initial_f = np.ones(days) * 0.01

#
res = minimize(objective, initial_f, method='L-BFGS-B', bounds=[(0, 2)]*days)

#
f_optimized = res.x

#
for t in range(1, days + 1):
    f = f_optimized[t-1]
    
    dw_dt[t] = f * a * b * W[t-1]**m
    W[t] = W[t-1] + dw_dt[t] * dt
    

#
def total_feeding_rate(t1, t2):
    return np.average(f_optimized[t1:t2])

print(total_feeding_rate(0, 114))

import pandas as pd
filename="fc.csv"
data = pd.read_csv(filename)
f1 = data['f'].to_numpy()

feed_1 = np.array([f_optimized[20],f_optimized[42],f_optimized[57],f_optimized[83],f_optimized[114]])

import matplotlib.pyplot as plt
plt.plot(f1[1:115], label='P', marker='o')
plt.plot(f_optimized, label='P_predict', marker='x')
plt.title('Comparison of P and P_predict')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.show()

Rmax = (17.477*W**(-0.39))/100
feed = []
for t in range(1, days + 1): 
    feed_d = f_optimized[t-1] * Rmax[t-1] * W[t-1]
    feed.append(feed_d)
    
feed2 = []
for t in range(1, days + 1): 
    feed_2 = f1[t] * Rmax[t-1] * W[t-1]
    feed2.append(feed_2)

FR = np.array(feed2) - np.array(feed)

plt.plot(FR, label='P', marker='o')
plt.title('Comparison of P and P_predict')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.show()






