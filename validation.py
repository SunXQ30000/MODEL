# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 16:25:28 2024

@author: Fighting！
"""

import numpy as np

# fixed parameters
#m = 0.54217624142879250470627994218375533819198608398438
m = 0.59662575
n = 5/6

DOcri= 0.6
DOmin= 0.3

UIAmax=1.452
UIAcri=0.1452

Topt = 30
Tmin = 15.8
Tmax = 46.8

a = 0.244491991

Ap = 0.9325833333
Af = 0.9817
Ac = 0.976725
  
kmin=0.008897945
j=0.037587236

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

Cp = 23.65655  # protein kj/g
Cf = 39.56715  # fats
Cc = 17.1667  # carbohydrates

Op = 1.89  # g O2/g protein
Of = 2.91  # g O2/g fat
Oc = 1.07  # g O2/g carbs

Pp = 0.1799

y = 0.8
Np = 1/6

#

def Rmax(W): 
    if W <=5:
        R = 0.1 
    elif 5<W<=25:
        R = 0.04
    elif 25<W<=100:
        R=0.03
    elif 100<W<=200:
        R=0.02
    elif 200<W<=300:
        R=0.015
    else:
        R=0.012
    return R

Fp_ = [0.535,0.531,0.516,0.487]
Ff_ = [0.1,0.11,0.125]
Fc = 0.08

Tact = 25
T = Tact

k = kmin * np.exp(j*(T-Tmin))

if T>Topt: 
    tau= np.exp( -4.6*((Topt-T)/(Topt-Tmin))**4 ) 
else:
    tau= np.exp( -4.6*((T-Topt)/(Tmax-Topt))**4 )

Tf = -0.0041204958*Tact**2 + 0.2570731016*Tact - 3.0379972353

def b(W):
    if W <= 10:
        Fp = Fp_[0]
        Ff = Ff_[0]
    elif 10 < W <= 100:
        Fp = Fp_[1]
        Ff = Ff_[0]
    elif 100 < W <= 300:
        Fp = Fp_[2]
        Ff = Ff_[1]
    else:
        Fp = Fp_[3]
        Ff = Ff_[2]
    
    return Fp*Ap + Ff*Af + Fc*Ac 

def F(W):
    if W <= 10:
        Fp = Fp_[0]
        Ff = Ff_[0]
    elif 10 < W <= 100:
        Fp = Fp_[1]
        Ff = Ff_[0]
    elif 100 < W <= 300:
        Fp = Fp_[2]
        Ff = Ff_[1]
    else:
        Fp = Fp_[3]
        Ff = Ff_[2]
    
    FP = 7.553958211*Fp**2 - 9.241008504*Fp + 3.721762897
    FF = -2.845528456*Ff + 1.287262873
    FC = 18.79699248*Fc**2 - 6.015037594*Fc + 1.360902256 
    return FP*FF*FC

"""
b = Fp*Ap + Ff*Af + Fc*Ac 

FP = 7.553958211*Fp**2 - 9.241008504*Fp + 3.721762897
FF = -2.845528456*Ff + 1.287262873
FC = 18.79699248*Fc**2 - 6.015037594*Fc + 1.360902256 
F = FP*FF*FC    
"""
# Initial conditions
W_initial = 11.6 # initial weight in grams
days = 300  # number of days
dt = 1  # time step in days
# Arrays to store results
weights = np.zeros(days + 1)
weights[0] = W_initial
Grow = np.zeros(days + 1)
Grow[0] = 0

for day in range(1, days + 1): 
    W = weights[day - 1]
    f = F(W) * Tf

    dw_dt = f * b(W) * tau * a * W ** m - k * W ** n
    W_new = W + dw_dt * dt

    weights[day] = W_new
    Grow[day] = dw_dt

from math import sqrt
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
W = weights
P = []
P_predict = []

rms1 = sqrt(mean_squared_error(P, P_predict))
print(rms1)
print(mean_absolute_error(P, P_predict))
print (r2_score(P, P_predict))

###############################
Pf = 3.1362* W ** 0.1652 /100

OCR = 0.003795798 * T - 0.01855184#0.000266364*T**2 - 0.009096602*T + 0.163924754 #湿重的OCR  kj/g

Fp = []
Ff = []
for i in range(len(W)):
    if W[i] <= 10:
        Fp0 = Fp_[0]
        Ff0 = Ff_[0]
    elif 10 < W[i] <= 100:
        Fp0 = Fp_[1]
        Ff0 = Ff_[0]
    elif 100 < W[i] <= 300:
        Fp0 = Fp_[2]
        Ff0 = Ff_[1]
    else:
        Fp0 = Fp_[3]
        Ff0 = Ff_[2]
    
    Fp.append(Fp0)
    Ff.append(Ff0)
Fp = np.array(Fp)
Ff = np.array(Ff)

SEC = (Fp*Cp)+(Ff*Cf)+(Fc*Cc) # KJ/g
#
Ep = (Fp*Cp)/SEC #protein
Ef = (Ff*Cf)/SEC #fats
Ec = (Fc*Cc)/SEC #carbohydrates

FL = ((1 - Ap) * Ep) + ((1 - Af) * Ef) + ((1 - Ac) * Ec)
BC = (0.494025609 * Ap * Ep) + (0.452951267 * Af * Ef) + (0.05 * Ac * Ec)
Cfi = Pp*Cp + Pf*Cf  # KJ g-1
e = 1 - FL - BC

## ====================
Qr = np.empty(days+1)
Qf = np.empty(days+1)
Qn = np.empty(days+1)
Qsda = np.empty(days+1)
Qg = np.empty(days+1)
Qs = np.empty(days+1)

Qg = Cfi * Grow
Qr = (OCR * weights **y + (Cfi - Np * Cp * Pp) * Grow) / (e - Np * Ep * Ap)
Qf= FL * Qr
Qn = Np * Cp * (Fp * Ap * (Qr / SEC) - Pp * Grow)
Qs = OCR * weights ** y
Qsda = BC * Qr

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
FCRt = Qr/(SEC*Grow)
feedt = Qr/SEC
def total_fcr(t1, t2):
    return np.average(FCRt[t1:t2])

print(total_fcr(1, 180))

feedt = Qr/SEC

output = {
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
