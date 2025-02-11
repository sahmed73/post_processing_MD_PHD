# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Fri Nov 29 11:32:19 2024
"""

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    y = np.zeros_like(x)
    
    y[(x >=0 ) & (x <= 1)] = 1
    y[(x > 1 ) & (x <= 2)] = -1

    return y

def Fourier(x,k):
    s=0
    amplitudes=[]
    frequencies=[]
    for i in range(k):
        n=2*i+1
        s+=(1/n)*np.sin(n*np.pi*x)
        frequencies.append(n/2)
        amplitudes.append(1/n)
        
    return s*4/np.pi, np.array(frequencies), np.array(amplitudes)
    
x = np.linspace(0, 2, 1000)
y = f(x)
Sn, frequencies, amplitudes = Fourier(x,10)

fig, ax = plt.subplots(1,2,dpi=350)
ax[0].plot(x,y,label='f(x)')
ax[0].plot(x,Sn,label='S$_n$')
ax[0].set_xlabel('x')
ax[0].set_ylabel('f(x) or S$_n$(x)')
ax[0].legend()
ax[1].scatter(frequencies,amplitudes,s=10,color='tab:red')
ax[1].set_xlabel('Frequency')
ax[1].set_ylabel('Amplitude')
plt.tight_layout()

