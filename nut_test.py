# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 12:30:24 2024

@author: Fighting！
"""

import numpy as np

# fixed parameters
#m = 0.54217624142879250470627994218375533819198608398438
m = 0.59662575

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

Fp_ = [0.4,0.5,0.6]
Ff_ = [0.3,0.2,0.1]
Fc = 0.08
Fp = Fp_[2]  #
Ff = Ff_[2]

Tact = 30
T = Tact

if T>Topt: 
    tau= np.exp( -4.6*((Topt-T)/(Topt-Tmin))**4 ) 
else:
    tau= np.exp( -4.6*((T-Topt)/(Tmax-Topt))**4 )
    
Tf = -0.0041204958*Tact**2 + 0.2570731016*Tact - 3.0379972353

b = Fp*Ap + Ff*Af + Fc*Ac 

FP = 7.553958211*Fp**2 - 9.241008504*Fp + 3.721762897
FF = -2.845528456*Ff + 1.287262873
FC = 18.79699248*Fc**2 - 6.015037594*Fc + 1.360902256 
F = FP*FF*FC    

# Initial conditions
W_initial = 1 # initial weight in grams
days = 660  # number of days
dt = 1  # time step in days
# Arrays to store results
weights = np.zeros(days + 1)
weights[0] = W_initial
Grow = np.zeros(days + 1)
Grow[0] = 0

for day in range(1, days + 1): 
    W = weights[day - 1]
    f = F*Tf       #feeding factor
    dw_dt = f * b * tau * a * W ** m
    W_new = W + dw_dt * dt

    weights[day] = W_new
    Grow[day] = dw_dt
"""
feed = []   
for t in range(1, days + 1):
    W = weights[t - 1]
    f = F*Tf
    feed_d = f * Rmax(W) * W
    feed.append(feed_d)
feed = np.array(feed)
FCR = []
for t in range(days): 
    FCRE = feed[t]/Grow[t+1]
    
    FCR.append(FCRE)
FCR = np.array(FCR)
def total_feeding_rate(t1, t2):
    return np.average(FCR[t1:t2])

print(total_feeding_rate(0, 360))
"""
################################  energy ############################
W = weights
Pf = 3.1362* W ** 0.1652 /100

OCR = 0.003795798 * T - 0.01855184

SEC = (Fp*Cp)+(Ff*Cf)+(Fc*Cc) # KJ/g
#
Ep = (Fp*Cp)/SEC #protein
Ef = (Ff*Cf)/SEC #fats
Ec = (Fc*Cc)/SEC #carbohydrates

FL = ((1 - Ap) * Ep) + ((1 - Af) * Ef) + ((1 - Ac) * Ec)
BC = (0.494025609 * Ap * Ep) + (0.452951267 * Af * Ef) + (0.05 * Ac * Ec)
Cfi = Pp*Cp + Pf*Cf  # KJ g-1
e = 1 - FL - BC

## ============= energy equation =============
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


#################################  density
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Parameters
initial_count = 10000
farm_area = 10*10*10  # m³
mortality_rate = 0.04
days_in = days

# Weight classes
weight_classes = ['<5', '5-10', '10-25', '25-100', '100-200', '200-300', '>300']
weight_limits = [5,10,25,100,200, 300]

# Arrays to track the number of sheep in each weight class and farming density
count_daily = np.zeros((days_in, 7))  # for each weight class
density_daily = np.zeros(days_in)

# Linear interpolation of sheep count due to mortality
survival_count = initial_count * (1 - mortality_rate)  # Final count after 6 months
count_interpolator = interp1d([0, days_in], [initial_count, survival_count], kind='linear')

# Simulate daily changes
for day in range(days_in):
    # Current sheep count based on linear interpolation
    current_count = count_interpolator(day)
    
    # Daily weight increase
    current_weight = W[day]
    
    # Classify sheep into weight classes based on current weight
    if current_weight < 5:
        count_daily[day, 0] = current_count
    elif current_weight < 10:
        count_daily[day, 1] = current_count
    elif current_weight < 25:
        count_daily[day, 2] = current_count
    elif current_weight < 100:
        count_daily[day, 3] = current_count
    elif current_weight < 200:
        count_daily[day, 4] = current_count
    elif current_weight < 300:
        count_daily[day, 5] = current_count
    else:
        count_daily[day, 6] = current_count
    
    # Calculate daily farming density (kg/m³)
    total_weight = current_count * current_weight
    density_daily[day] = total_weight / farm_area

# Plotting sheep count changes over time in each weight class
plt.rcParams['font.family'] = 'Times New Roman'

plt.figure(figsize=(10, 6))
for i in range(4):
    plt.plot(count_daily[:, i], label=weight_classes[i])
plt.xlabel('Days')
plt.ylabel('Sheep Count')
plt.tick_params(axis='both', direction='in')
plt.legend()
plt.show()

# Plotting daily farming density
days = np.arange(len(density_daily))  # Assuming you have the number of days as the x-axis
plt.figure(figsize=(10, 6))
plt.fill_between(days, density_daily, where=density_daily > 200, color='lightcoral', label='Density > 30 (kg/m³)')
plt.fill_between(days, density_daily, where=density_daily <= 200, color='lightblue', label='Density <= 30 (kg/m³)')
plt.plot(density_daily, label='Farming Density (kg/m³)', color='r')
plt.xlabel('Days')
plt.ylabel('Density (kg/m³)')
plt.tick_params(axis='both', direction='in')
plt.ylim(100, 400)
plt.show()



