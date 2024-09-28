# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 12:32:07 2024

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
s_ = [0.8, 0.9, 1.0]
s = s_[2]  #

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
FC = 18.79699248*Fc**2 - 6.015037594*Fc + 1.360902256  

Tact = 30

T = Tact

if T>Topt: 
    tau= np.exp( -4.6*((Topt-T)/(Topt-Tmin))**4 ) 
else:
    tau= np.exp( -4.6*((T-Topt)/(Tmax-Topt))**4 )

Tf = -0.0041204958*T**2 + 0.2570731016*T - 3.0379972353

# Initial conditions
W_initial = 1 # initial weight in grams
days = 660  # number of days
dt = 1  # time step in days
# Arrays to store results
weights = np.zeros(days + 1)
weights[0] = W_initial
Grow = np.zeros(days + 1)
Grow[0] = 0

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
    
    return FP*FF*FC

for day in range(1, days + 1): 
    W = weights[day - 1]
    #current_r = r[(day-1) // switch_interval]
    f = s * F(W) * Tf       
    dw_dt = f * b(W) * tau * a * W ** m
    W_new = W + dw_dt * dt

    weights[day] = W_new
    Grow[day] = dw_dt
"""
feed = []   
for t in range(1, days + 1):
    W = weights[t - 1]
    f = s * F(W) * Tf
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
################################
W = weights
Pf = 3.1362* W ** 0.1652 /100

OCR = 0.003795798 * T - 0.01855184

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

## =============  =============
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





