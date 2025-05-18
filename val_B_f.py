# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize

# Parameters
m = 0.6121

Fp_ = [0.49,0.49,0.49,0.49]  #[0.535,0.531,0.516,0.487]
Ff_ = [0.09,0.09,0.09] #  [0.1,0.11,0.125]
Fc = 0.08

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

a = 0.2441

Topt = 30
Tmin = 15.8
Tmax = 46.8
Tact = 26

Tf = -0.0041204958*Tact**2 + 0.2570731016*Tact - 3.0379972353

dt = 1  # time step in days

W_observed = np.array([156.3,198.5,274.7,302.4,563.4,792.3,1026.6,1324.9])  # 
t_observed = np.array([0,31,63,94,126,158,189,221])  #

W0 = W_observed[0]  # 
days = t_observed[-1]  # 
W = np.zeros(days + 1)  # 
W[0] = W0
dw_dt = np.zeros(days + 1)  # 
feeding_rates = np.zeros(days + 1)  # 

# 
def objective(f_values):
    W_temp = np.zeros_like(W)
    W_temp[0] = W0
    for t in range(1, len(W_temp)):
        f = f_values[t-1]  # 
        dw_dt = f * a * b(W_temp[t-1]) * W_temp[t-1]**m
        W_temp[t] = W_temp[t-1] + dw_dt * dt
    return np.sum((W_temp[t_observed] - W_observed)**2)  # 

# 
initial_f =  np.full(days, 0.5)

# 
res = minimize(objective, initial_f, method='L-BFGS-B', bounds=[(0, 2)]*days)

# 
f_optimized = res.x

# 
for t in range(1, days + 1):
    f = f_optimized[t-1]
    dw_dt[t] = f * a * b(W[t-1]) * W[t-1]**m
    
    W[t] = W[t-1] + dw_dt[t] * dt


def total_feeding_rate(t1, t2):
    return np.average(f_optimized[t1:t2])

print(total_feeding_rate(0, 220))

P1 = [0.410245523,0.558324847,0.175931185,1.460935144,0.897890328,0.854002553,0.864827751]
P_predict = [f_optimized[30],f_optimized[62],f_optimized[93],f_optimized[125],f_optimized[157],f_optimized[188],f_optimized[220]]
f_predict = np.array(P_predict)

import matplotlib.pyplot as plt
plt.plot(P1, label='P', marker='o')
plt.plot(P_predict, label='P_predict', marker='x')
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
    
import pandas as pd
filename="fb.csv"
data = pd.read_csv(filename)
P = data['f'].to_numpy()
feed2 = []
for t in range(1, days + 1): 
    feed_2 = P[t-1] * Rmax[t-1] * W[t-1]
    feed2.append(feed_2)

FR = np.array(feed2) - np.array(feed)

plt.plot(FR, label='P', marker='o')
plt.title('Comparison of P and P_predict')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.show()

################################################################
Grow = dw_dt
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

Pf = (3.1362* W ** 0.1652) /100

OCR = 0.003795798 * Tact - 0.01855184

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
Qs = OCR * W ** y
Qsda = BC * Qr

EN = np.empty(len(W))
Efae = np.empty(len(W))
DO2 = np.empty(len(W))
DO2fe = np.empty(len(W))
ER = np.empty(len(W))

MP = Fp * Ap * Qr / SEC - Pp * Grow
ML = Ff * Af * Qr / SEC - Pf * Grow
CAR = Fc * Ac * Qr / SEC

DO2fr = Fp*Op + Ff*Of + Fc*Oc

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













