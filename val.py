# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 19:11:40 2024

@author: Fighting！
"""

import numpy as np

# fixed parameters
m = 0.61
n = 5/6

DOcri= 0.6
DOmin= 0.3

UIAmax=1.452
UIAcri=0.1452

Topt = 30
Tmin = 15.8
Tmax = 46.8

a = 0.24

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

#非固定
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

Fp = 0.465  #0.45-0.48  0.08-0.1
Ff = 0.09
Fc = 0.08

Tact = [24.3,24.3,24.5,24.6,24.7,25.5,25.8,26.1,26.7,27.2,28.5,28.1]
T = Tact

k=[]*len(Tact)
for i in range(len(Tact)): 
    k.append(kmin * np.exp(j*(Tact[i]-Tmin)))

Tf=[]*len(Tact)
for i in range(len(Tact)): 
    if Tact[i]>Topt: 
        tau= np.exp( -4.6*((Topt-Tact[i])/(Topt-Tmin))**4 ) 
    else:
        tau= np.exp( -4.6*((Tact[i]-Topt)/(Tmax-Topt))**4 )
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
Grow[0] = 0

switch_interval = 15
for day in range(1, days + 1): 
    W = weights[day - 1]
    current_T = T[(day-1) // switch_interval]
    current_Tf = Tf[(day-1) // switch_interval]
    #* current_Tf* F
    #current_r = r[(day-1) // switch_interval]
    f = current_Tf
    #f= r(W)/Rmax(W)
    dw_dt = b * a* current_T * W ** m
    W_new = W + dw_dt * dt

    weights[day] = W_new
    Grow[day] = dw_dt
    
"""
feed = []   
for t in range(1, days + 1): 
    W = weights[t-1]
    current_Tf = Tf[(day-1) // switch_interval]
    #* current_Tf* F
    #current_r = r[(day-1) // switch_interval]
    f = F * current_Tf
    feed_d = f * Rmax(W) * W
    feed.append(feed_d)
    #* current_Tf* F 
FCR = []
for t in range(days): 
    FCRE = feed[t]/Grow[t+1]
    
    FCR.append(FCRE)
"""
from math import sqrt
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
W = weights
P = [10.9,14.9,21.3,29.9,50.4,70.6,100.3,141.2,201.4,286.4,370.4,460,571.2]
P_predict = [W[0],W[14],W[29],W[44],W[59],W[74],W[89],W[104],W[119],W[134],W[149],W[164],W[179]]

rms1 = sqrt(mean_squared_error(P, P_predict))
print(rms1)
print(mean_absolute_error(P, P_predict))
print (r2_score(P, P_predict)) 

################################  energy
Pf = 3.1362* W ** 0.1652 /100

OCR1=[]*len(Tact)
for i in range(len(Tact)): 
    OCR_ = 0.003795798 * Tact[i] - 0.01855184      #0.000266364*Tact[i]**2 - 0.009096602*Tact[i] + 0.163924754 #湿重的OCR  kj/g
    OCR1.append(OCR_)
    
OCR = np.empty(days+1)
for t in range(days): 
    OCR[t] = OCR1[(day-1) // switch_interval]

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

Qf = np.empty(days+1)
Qn = np.empty(days+1)
Qsda = np.empty(days+1)
Qg = np.empty(days+1)
Qs = np.empty(days+1)
Qr = np.empty(days+1)

Qg = Cfi * Grow
Qr = (OCR * weights **y + (Cfi - Np * Cp * Pp) * Grow) / (e - Np * Ep * Ap)
Qf= FL * Qr
Qn = Np * Cp * (Fp * Ap * (Qr / SEC) - Pp * Grow)
#Qs = OCR * weights ** y
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

print(total_fcr(1, 181))

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

#################################  density
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Parameters
initial_count = 10000
farm_area = 40*1000  # m³
mortality_rate = 0.0262
days_in = days + 1

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
for i in range(7):
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
plt.fill_between(days, density_daily, where=density_daily <= 300, color='lightblue', label='Density <= 30 (kg/m³)')
plt.plot(density_daily, label='Farming Density (kg/m³)', color='r')
plt.xlabel('Days')
plt.ylabel('Density (kg/m³)')
plt.tick_params(axis='both', direction='in')
plt.ylim(0, 150)
plt.show()


