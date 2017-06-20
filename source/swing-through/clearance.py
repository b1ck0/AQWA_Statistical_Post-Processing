from math import pi, sqrt, log, gamma, cos, sin

from matplotlib.backends.backend_pdf import PdfPages

from scipy.stats import norm, rayleigh
from scipy.signal import welch

import numpy as np

import matplotlib.pyplot as plt

class Globs:
    FS = 1.0 #sampling frequency
    N = 256 #overlaping intervals
    t = 10800 #seconds
    L = 270 #chain length in [m]
         
def printExtremeValue(array,FS,N):
    xf, yf = welch(array,fs=FS,window='hann',nperseg=N,return_onesided=True)
    
    return Rmax(xf,yf)+np.mean(array)
                
def Rmax(xf,yf):
    t = Globs.t
    yf2 = yf * xf**2 #ordinates of the second order spectrum
    M0 = np.trapz(yf,xf) #integrating the area under the yf curve
    M2 = np.trapz(yf2,xf) #integrating the area under the yf2 curve
    std = np.sqrt(M0) #obtaining the standard deviation of the input signal = square root of the variance
    Vx = np.sqrt(M2/M0) #average zero-upcrossing frequency of the input signal
    p = 0.37 #probability of not exceedence of the Most Probable Extreme
    Max = std*np.sqrt(2.0*np.log((Vx*t)/(np.log(1.0/p))))
    
    return Max
       
def readCSV(filename):
    values = []
    
    with open(filename) as f:
        lines = f.readlines()
    
    for ID in range(7,len(lines)):
        strList = lines[ID].split(",")
        #storing only the whole and half seconds
        if ((float(strList[0])).is_integer()): #or ((float(strList[0])-0.5).is_integer()):
            values.append(float(strList[1]))
    
    return values
    
def getSign(a):
    if a > 0:
        return 1
    else: 
        return -1

#execution of the main module

P1X = []
P1Y = []
P2X = []
P2Y = []
D1 = []
D2 = []

p1x = []
p1y = []
p2x = []
p2y = []

P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA046.CSV")))
P1X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA046.CSV")))

P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA047.CSV")))
P1Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA047.CSV")))

P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA048.CSV")))
P2X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA048.CSV")))

P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA049.CSV")))
P2Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA049.CSV")))

p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA050.CSV")))
p1x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA050.CSV")))

p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA051.CSV")))
p1y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA051.CSV")))

p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA052.CSV")))
p2x.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA052.CSV")))

p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA053.CSV")))
p2y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA053.CSV")))

#----------------------------------------------------------------------reading the result files 

D1 = []
D2 = []

maxD1 = []
maxD2 = []

for i in range(0,len(P1X)):
    d1 = []
    d2 = []
    for j in range(0,len(P1X[i])):
        d1.append(getSign(P1X[i][j] - p1x[i][j])*sqrt((P1X[i][j] - p1x[i][j])**2 + (P1Y[i][j] - p1y[i][j])**2))
        d2.append(getSign(P2X[i][j] - p2x[i][j])*sqrt((P2X[i][j] - p2x[i][j])**2 + (P2Y[i][j] - p2y[i][j])**2))
    D1.append(d1)
    D2.append(d2)

print len(D1)
print len(D2)

for i in D1:
    maxD1.append(round(min(i),2))
    
for i in D2:
    maxD2.append(round(max(i),2))

print maxD1
print maxD2