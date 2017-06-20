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
    array2 = []
    
    for i in array:
        array2.append(270.0-i)
            
    xf, yf = welch(array2,fs=FS,window='hann',nperseg=N,return_onesided=True)
    
    return Rmax(xf,yf)+np.mean(array2)
                
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

STB = []
PS = []

STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA009.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA009.CSV")

PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA012.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA012.CSV")

#----------------------------------------------------------------------reading the result files 

resSTB = []
resPS = []

maxSTB = []
maxPS = []

for i in range(0,8,1):
    resSTB.append(readCSV(STB[i]))
    resPS.append(readCSV(PS[i]))

limit = 1386.0    
            
for i in range(0,8,1):
    maxSTB.append(round(270-printExtremeValue(resSTB[i],Globs.FS,Globs.N),1))
    maxPS.append(round(270-printExtremeValue(resPS[i],Globs.FS,Globs.N),1))
    
print maxSTB
print maxPS