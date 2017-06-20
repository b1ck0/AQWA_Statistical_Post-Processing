from math import pi, sqrt, log, gamma, cos, sin, radians

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
        #if ((float(strList[0])).is_integer()): #and ((float(strList[0])) > 0): #or ((float(strList[0])-0.5).is_integer()):
        values.append(float(strList[1]))
    
    return values

def FSTYAW_1(yaw, Rx, Ry, sign): #transofrmation of axis sing = 1 --> FRA to LSA; sign == -1 --> LSA to FRA
    roll = 0.0
    pitch = 0.0
    
    c = cos(radians(yaw))
    s = sin(radians(yaw))
    
    if sign == 1:
        roll = round(Rx*c + Ry*s,1)
        pitch = round(-Rx*s + Ry*c,1)
    else:
        roll = round(Rx*c - Ry*s,1)
        pitch = round(Rx*s + Ry*c,1)
        
    return roll, pitch             
       
def FRAtoLSA(Rx,Ry,Rz):
    roll, pitch = [], []
    for i in range(0,len(Rx)):
        r, p = FSTYAW_1(Rz[i],Rx[i],Ry[i],1)
        roll.append(r)
        pitch.append(p)
        
    return roll,pitch
#execution of the main module

X = []
Y = []
Z = []
YAW = []
F = []

X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA007.CSV"))) #7 or 10
#X.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA010.CSV"))) #7 or 10

Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA008.CSV"))) #8 or 11
#Y.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA011.CSV"))) #8 or 11

Z.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA009.CSV"))) #9 or 12
#Z.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA012.CSV"))) #9 or 12

YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5-0117 - Husky/results-swing through/Tp11_R2/AQWA006.CSV")))

#----------------------------------------------------------------------reading the result files 

xLoc = []
yLoc = []

xLocT, yLocT = FRAtoLSA(X[0],Y[0],YAW[0])

for i in range(0,len(xLocT)):
    xLoc.append(abs(xLocT[i]))
    yLoc.append(abs(yLocT[i]))
    
#xLoc.append(xLocT)
#yLoc.append(yLocT)
   
for i in range(0,len(X[0])):
    f = round(sqrt(X[0][i]**2+Y[0][i]**2+Z[0][i]**2),1)
    F.append(f)
   
print round(max(xLoc)/1000.0,1)
print round(max(yLoc)/1000.0,1)
print round(max(Z[0])/1000.0,1)
print round(max(F)/1000.0,1)