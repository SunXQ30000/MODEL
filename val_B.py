# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 16:25:15 2024

@author: Fighting！
"""

import numpy as np

# fixed parameters
m = 0.61

DOcri= 0.6
DOmin= 0.3

UIAmax=1.452
UIAcri=0.1452
"""
Topt = 30
Tmin = 15.8
Tmax = 46.8
"""
Topt = 30
Tmin = 14
Tmax = 40

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

#非固定

Fp_ = [0.49,0.49,0.49,0.49] ##[0.535,0.531,0.516,0.487]
Ff_ = [0.09,0.09,0.09]  #[0.1,0.11,0.125]
Fc = 0.08

Tact = 26.5  #22.3-33.726.5
T = Tact

#k = kmin * np.exp(j*(T-Tmin))
if T<Topt: 
    tau= np.exp( -4.6*((Topt-T)/(Topt-Tmin))**4 ) 
else:
    tau= np.exp( -4.6*((T-Topt)/(Tmax-Topt))**4 )

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
    
    return Fp + Ff + Fc 

# Initial conditions
W_initial = 80 # initial weight in grams
days = 115  # number of days
dt = 1  # time step in days
# Arrays to store results
weights = np.zeros(days + 1)
weights[0] = W_initial
Grow = np.zeros(days + 1)
Grow[0] = 0
feed = np.zeros(days + 1)
feed[0] = 0

import pandas as pd
filename="fc.csv"
data = pd.read_csv(filename)
f1 = data['f2'].to_numpy()

for day in range(1, days + 1): 
    W = weights[day - 1]
    #f = f1[day]
    
    f = 0.763945743*tau
    
    dw_dt = f * b(W) * a * W ** m
    W_new = W + dw_dt * dt

    weights[day] = W_new
    Grow[day] = dw_dt
    feed[day] = f

from math import sqrt
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
W = weights
P = np.array([80,130.89736,174.4584,204.76174,334.60311,431])
P_predict = np.array([W[0],W[20],W[42],W[57],W[83],W[114]])
feed_1 = [feed[0],feed[20],feed[42],feed[57],feed[83],feed[115]]
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

OCR = 0.003795798 * T - 0.01855184

Fp = 0.49
Ff = 0.09

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







