# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 19:11:40 2024

@author: Fighting！
"""

import numpy as np

# fixed parameters
m = 0.61

DOcri= 0.6
DOmin= 0.3

UIAmax=1.452
UIAcri=0.1452

Topt = 30
Tmin = 15.8
Tmax = 46.8

a = 0.24

DO = 7
UIA = 0.001

if DO>DOcri: 
    sigma=1
elif DO>DOmin: 
    sigma= (DO-DOmin)/(DOcri-DOmin)
else:
    sigma=0
    
if UIA<UIAcri: 
    v=1
elif UIA<UIAmax: 
    v= (UIAmax-UIA)/(UIAmax-UIAcri)
else: 
    v=0

#

Fp = 0.47  #0.45-0.48  0.08-0.1
Ff = 0.09
Fc = 0.08

Tact = [24.3,24.3,24.5,24.6,24.7,25.5,25.8,26.1,26.7,27.2,28.5,28.1]
T = Tact 
    
Tf=[]*len(Tact)
for i in range(len(Tact)): 
    if Tact[i]<Topt: 
        T = Tact[i]
        tau= np.exp( -4.6*((Topt-T)/(Topt-Tmin))**4 ) 
    else:
        tau= np.exp( -4.6*((T-Topt)/(Tmax-Topt))**4 )
    Tf.append(tau)
    
b = Fp + Ff + Fc

# Initial conditions
W_initial = 10.9 # initial weight in grams
days = 180  # number of days
dt = 1  # time step in days
# Arrays to store results
weights = np.zeros(days + 1)
weights[0] = W_initial
Grow = np.zeros(days + 1)
feed = np.zeros(days + 1)
Grow[0] = 0
feed[0]= 0

#
f1 = [0.363,0.452329,0.502585,0.918745,0.696974,0.829973,0.924013,1.09984,1.25015,1.02552,0.945774,1.02823]
#
#f1 = [0.362343576,0.471537614,0.512852374,0.929021947,0.71285449,0.84904828,0.946874566,1.125964759,1.281563193,1.056315527,0.976542902,1.0619748]

switch_interval = 15
for day in range(1, days + 1): 
    W = weights[day - 1]

    f = 0.9288086 * Tf[(day-1) // switch_interval]
    
    #f = f1[(day-1) // switch_interval]
    
    dw_dt = f * b * a * W ** m
    W_new = W + dw_dt * dt
    
    feed[day] = f
    weights[day] = W_new
    Grow[day] = dw_dt
    
from math import sqrt
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
W = weights
P = np.array([10.9,14.9,21.3,29.9,50.4,70.6,100.3,141.2,201.4,286.4,370.4,460,571.2])
P_predict = np.array([W[0],W[14],W[29],W[44],W[59],W[74],W[89],W[104],W[119],W[134],W[149],W[164],W[179]])
feed_1 = [feed[14],feed[29],feed[44],feed[59],feed[74],feed[89],feed[104],feed[119],feed[134],feed[149],feed[164],feed[179]]
feed_1 = np.array(feed_1)
w_predict = np.array(P_predict)

rms1 = sqrt(mean_squared_error(P, P_predict))
print(rms1)
mape = np.mean(np.abs((P - P_predict) / P)) * 100
print(mape)
print(mean_absolute_error(P, P_predict))
print (r2_score(P, P_predict)) 

import matplotlib.pyplot as plt
plt.plot(P, label='P', marker='o')
plt.plot(P_predict, label='P_predict', marker='x')
plt.title('Comparison of P and P_predict')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.show()

################################################################
Cp = 23.65655  # protein kj/g
Cf = 39.56715  # fats
Cc = 17.1667  # carbohydrates

Op = 1.89  # g O2/g protein
Of = 2.91  # g O2/g fat
Oc = 1.07  # g O2/g carbs

Pp = 0.1799

y = 0.8
Np = 1/6

Ap = 0.9325833333
Af = 0.9817
Ac = 0.976725

W = weights
Pf = (3.1362* W ** 0.1652) /100

OCR=[]*len(W)
for i in range(len(W)):
    OCR_ = 0.003795798 * Tact[(i-1) // switch_interval] - 0.01855184
    OCR.append(OCR_)

SEC = (Fp*Cp)+(Ff*Cf)+(Fc*Cc) # KJ/g
#
Ep = (Fp*Cp)/SEC #protein
Ef = (Ff*Cf)/SEC #fats
Ec = (Fc*Cc)/SEC #carbohydrates

FL = ((1 - Ap) * Ep) + ((1 - Af) * Ef) + ((1 - Ac) * Ec)
#BC = (0.494025609 * Ap * Ep) + (0.452951267 * Af * Ef) + (0.05 * Ac * Ec)
BC = 0.3 * Ap * Ep + 0.05 * Af * Ef + 0.05 * Ac * Ec
Cfi = Pp*Cp + Pf*Cf  # KJ g-1
e = 1 - FL - BC

## ============= energy equation =============
Qr = np.empty(days+1)
Qf = np.empty(days+1)
Qn = np.empty(days+1)
Qsda = np.empty(days+1)
Qg = np.empty(days+1)
Qs = np.empty(days+1)

Rmax = (17.477*W**(-0.39))/100
food = Rmax*W*feed

Qg = Cfi * Grow
#Qr = (OCR * weights **y + (Cfi - Np * Cp * Pp) * Grow) / (e - Np * Ep * Ap)
Qr = food * SEC
Qf= FL * Qr
Qn = Np * Cp * (Fp * Ap * (Qr / SEC) - Pp * Grow)
Qs = OCR * weights ** y
Qsda = BC * Qr

DO2fr = Fp*Op + Ff*Of + Fc*Oc

EN = np.empty(len(W))
Efae = np.empty(len(W))
DO2 = np.empty(len(W))
DO2fe = np.empty(len(W))
ER = np.empty(len(W))

MP = Fp * Ap * Qr / SEC - Pp * Grow
ML = Ff * Af * Qr / SEC - Pf * Grow
CAR = Fc * Ac * Qr / SEC

EN = MP * Np
Efae = FL * Qr / SEC
DO2 = MP * Op + ML * Of + CAR * Oc
DO2fe = (Fp*(1-Ap)*Op + Ff*(1-Af)*Of + Fc*(1-Ac) * Oc) * Qr / SEC
ER = Qg / Qr
FCRt = food/Grow
feedt = Qr/SEC
def total_fcr(t1, t2):
    return np.average(FCRt[t1:t2])

print(total_fcr(1, 180))

feedt = Qr/SEC

QQQQQ = {
    'Qr': Qr,
    'Qg': Qg,
    'Qf': Qf,
    'Qn': Qn,
    'Qs': Qs,
    'Qsda': Qsda
}

aoutput = {
    'weight': W,
    'FCRt': FCRt,
    'feedt': feedt,
    'ER': ER,
    'EN': EN,
    'Efae': Efae,
    'DO2': DO2,
    'DO2fe': DO2fe
}

#################################
num_chicks = 10000  # total number of chicks
survival_rate = 0.9684  # survival rate after 200 days
pool_volume = 5 * 5 * 1.5  # volume of the pond in cubic meters

# Weight categories in grams
categories = {
    "0-10": (0, 10),
    "10-100": (10, 100),
    "100-300": (100, 300),
    "300+": (300, 4000),
}

# Generate weight growth data
np.random.seed(42)  # for reproducibility
daily_weights = weights

# Simulate daily survival
daily_survival = np.linspace(1, survival_rate, days)
num_alive = num_chicks * daily_survival

# Categorize chicks by weight daily
category_counts = {cat: [] for cat in categories}

total_biomass = []  # total biomass (in grams) of surviving chicks
daily_density = []  # density in kg/m^3
for i in range(days):
    weights_today = np.random.normal(daily_weights[i], 50, int(num_alive[i]))
    total_biomass_today = np.sum(weights_today)  # total biomass in grams
    total_biomass.append(total_biomass_today)

    # Calculate density in kg/m^3 (convert grams to kg)
    density_today = (total_biomass_today / 1000) / pool_volume  # in kg/m^3
    daily_density.append(density_today)
    
    for cat, (low, high) in categories.items():
        count_in_cat = np.sum((weights_today >= low) & (weights_today < high))
        category_counts[cat].append(count_in_cat)


import matplotlib.pyplot as plt

plt.figure(figsize=(12, 6))
# Plot weight categories over time
for cat, counts in category_counts.items():
    plt.plot(counts, label=cat)

plt.xlabel('Days')
plt.ylabel('Number of Chicks')
plt.title('Number of Chicks in Different Weight Categories Over Time')
plt.legend()
plt.grid(True)
plt.show()



