# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:00:46 2025

@author: 97208
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from SALib.sample import saltelli
from SALib.analyze import sobol

# fixed parameters
#m = 0.54217624142879250470627994218375533819198608398438
m = 0.6121

DOcri= 0.6
DOmin= 0.3

UIAmax=1.452
UIAcri=0.1452

a = 0.2441

def model(params):
    T, Fp, Ff, Fc, s, DO, UIA = params        #
    
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
    
    Tf = -0.0041204958*T**2 + 0.2570731016*T - 3.0379972353
    
    b = Fp + Ff + Fc

    FP = 7.553958211*Fp**2 - 9.241008504*Fp + 3.721762897
    FF = -2.845528456*Ff + 1.287262873
    FC = 18.79699248*Fc**2 - 6.015037594*Fc + 1.360902256 
    F = FP*FF*FC

    W_initial = 10 # initial weight in grams
    days = 365  # number of days
    dt = 1  # time step in days

    weights = np.zeros(days + 1)
    weights[0] = W_initial
    Grow = np.zeros(days + 1)
    feed = np.zeros(days + 1)
    Grow[0] = 0
    feed[0]= 0

    switch_interval = 15
    for day in range(1, days + 1): 
        W = weights[day - 1]

        f = s *F * Tf *sigma *v
        
        dw_dt = f * b * a * W ** m
        W_new = W + dw_dt * dt
    
        feed[day] = f
        weights[day] = W_new
        Grow[day] = dw_dt
        
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

    # print(total_fcr(1, 180))

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
    
    return sum(feedt[1:]) #### FCRt[-1]  / len(FCRt[1:])         sum(numbers[1:]) / len(numbers[1:])
    
#非固定
problem = {
    'num_vars': 7,
    'names': ['T', 'Fp', 'Ff', 'Fc', 's', 'DO', 'UIA'],   #
    'bounds': [
        [16, 35],    # T范围
        [0.45, 0.55],     # DO范围0.3,0.45
        [0.07, 0.13],     # UIA范围（保证FUIA非负）
        [0.08, 0.16],
        [0, 1],
        [0.3, 0.6],           #  [5, 8]
        [0.1452, 1.452]  # N范围  [0.001, 0.01]
    ]
}


# 生成Sobol样本
param_values = saltelli.sample(problem, 1024)

Y = np.array([model(x) for x in param_values])

# 计算灵敏度指数
Si = sobol.analyze(problem, Y)

# 输出结果
print(Si['S1'])  # 主效应
print(Si['ST'])  # 总效应
print(Si['S2'])  # 交互效应

interaction_matrix = Si['S2']
interaction_matrix = np.array(interaction_matrix).reshape((7, 7))

# 使用 seaborn 绘制热力图
plt.figure(figsize=(10, 8))
sns.heatmap(interaction_matrix, annot=True, fmt=".4f", cmap="YlGnBu", xticklabels=problem['names'], yticklabels=problem['names'])

# 设置标题和显示热力图
plt.title('Parameter Interaction Effects')
plt.show()


