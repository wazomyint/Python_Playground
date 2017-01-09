# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:39:33 2016
NMR peak creator
@author: wazomyint
"""

import tkinter as tk
import numpy as np
import math
import matplotlib.pyplot as plt
from tkinter.filedialog import asksaveasfilename

def callback():
    freq1 = float(E1.get())*2*np.pi
    freq2 = float(E2.get())*2*np.pi
    R2 = float(E6.get())
    acq = float(E3.get())
    dw = float(E4.get())*1e-6
    err = float(E5.get())*0.01
    temp = acq/dw
    Jcoup= float(E6.get())*2*np.pi
    xpts = np.linspace(0, acq, temp)
    ypts = []
    for x in xpts:
         ypts.append((math.cos((freq1+Jcoup/2)*x)+1j*math.sin((freq1+Jcoup/2)*x)+(math.cos((freq1-Jcoup/2)*x)+1j*math.sin((freq1-Jcoup/2)*x))+ \
         (math.cos((freq2+Jcoup/2)*x)+1j*math.sin((freq2+Jcoup/2)*x)+(math.cos((freq2-Jcoup/2)*x)+1j*math.sin((freq2-Jcoup/2)*x))))*math.exp(-R2*x))
    ypts += err * np.random.randn(len(ypts))
    fidfig = plt.figure()
    p = fidfig.add_subplot(1, 1, 1)
    line_test, = p.plot(xpts, ypts)
    ftfig = plt.figure()
    ft = ftfig.add_subplot(1,1,1)
    spec = np.fft.fft(ypts)
    freq = np.fft.fftfreq(len(ypts), dw)
    F_shifted = np.fft.fftshift(spec)
    freq_shifted = np.fft.fftshift(freq)
    line_test, = ft.plot(freq_shifted, np.abs(F_shifted))   
    savename = asksaveasfilename(title = "Save Peak in EPS")
    plt.savefig(savename, format = 'eps')
    top.quit()
    top.destroy()
    
def _quit():
    top.quit()
    top.destroy()
    
# Main Window
top = tk.Tk()

# Freq 1
L1 = tk.Label(top, text = "Frequency 1 (Hz)")
L1.grid(row=0, column=0)
E1 = tk.Entry(top, bd=3)
E1.grid(row=0, column=1)

# Freq 2
L2 = tk.Label(top, text = "Frequency 2 (Hz)")
L2.grid(row=1, column=0)
E2 = tk.Entry(top, bd=3)
E2.grid(row=1, column=1)

# Total Acquisition
L3 = tk.Label(top, text = "Total Acquisition (s)")
L3.grid(row=2, column=0)
E3 = tk.Entry(top, bd=3)
E3.grid(row=2, column=1)

# Dwell time
L4 = tk.Label(top, text = "Dwell Time (\u03bcs)")
L4.grid(row=3, column=0)
E4 = tk.Entry(top, bd=3)
E4.grid(row=3, column=1)

# Noise
L5 = tk.Label(top, text = "Error as %")
L5.grid(row=4, column=0)
E5 = tk.Entry(top, bd=3)
E5.grid(row=4, column=1)

# Relaxation Rate
L6 = tk.Label(top, text = "Relaxation Rate")
L6.grid(row=5, column=0)
E6 = tk.Entry(top, bd=3)
E6.grid(row=5, column=1)

# J coupling
L7 = tk.Label(top, text = "J Coupling")
L7.grid(row=6, column=0)
E7 = tk.Entry(top, bd=3)
E7.grid(row=6, column=1)


B1 = tk.Button(top, text = "Make Graph", command =callback)
B1.grid(row = 7, column  = 0)
B2 = tk.Button(top, text = "Quit", command = _quit)
B2.grid(row = 7, column = 1)
top.mainloop()
