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

#execution of the main module

X = []
Y = []

X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp4/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp5/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp6/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp7/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp8/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp9/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp10/AQWA044.CSV")))
X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp11/AQWA044.CSV")))

Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp4/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp5/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp6/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp7/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp8/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp9/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp10/AQWA045.CSV")))
Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp11/AQWA045.CSV")))

#----------------------------------------------------------------------reading the result files 

TOTAL = []

maxX = []
maxY = []
maxT = []

for i in range(0,len(X)):
    arr = []
    for j in range(0,len(X[i])):
        t = round(sqrt(X[i][j]**2+Y[i][j]**2),2)
        arr.append(t)
    TOTAL.append(arr)

for i in range(0,len(X)):
    maxX.append(round(max(X[i])/1000.0,1))
    maxY.append(round(max(Y[i])/1000.0,1))
    maxT.append(round(max(TOTAL[i])/1000.0,1))
    
print maxX
print maxY
print maxT