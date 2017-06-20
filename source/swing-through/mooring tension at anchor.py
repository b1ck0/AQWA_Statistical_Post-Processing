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
        if ((float(strList[0])).is_integer()) and ((float(strList[0])) > 1200): #or ((float(strList[0])-0.5).is_integer()):
            values.append(float(strList[1]))
    
    return values
       
#execution of the main module

STB = []
PS = []

STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA054.CSV")
STB.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA054.CSV")

PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp4/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp5/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp6/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp7/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp8/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp9/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp10/AQWA055.CSV")
PS.append("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_swingthrough/Tp11/AQWA055.CSV")

#----------------------------------------------------------------------reading the result files 

resSTB = []
resPS = []

maxSTB = []
maxPS = []

for i in range(0,8,1):
    resSTB.append(readCSV(STB[i]))
    resPS.append(readCSV(PS[i]))

limit = 676.0    
            
for i in range(0,8,1):
    maxSTB.append(round(printExtremeValue(resSTB[i],Globs.FS,Globs.N)/1000.0/limit,2))
    maxPS.append(round(printExtremeValue(resPS[i],Globs.FS,Globs.N)/1000.0/limit,2))
    
print maxSTB
print maxPS