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
        if ((float(strList[0])).is_integer()): #or ((float(strList[0])-0.5).is_integer()):
            values.append(float(strList[1]))
    
    return values

def FSTYAW_1(yaw, Rx, Ry, sign): #transofrmation of axis sing = 1 --> FRA to LSA; sign == -1 --> LSA to FRA
    roll = 0.0
    pitch = 0.0
    
    c = cos(radians(yaw))
    s = sin(radians(yaw))
    
    if sign == 1:
        roll = Rx*c + Ry*s
        pitch = -Rx*s + Ry*c
    else:
        roll = Rx*c - Ry*s
        pitch = Rx*s + Ry*c
        
    return roll, pitch             
       
def FRAtoLSA(Rx,Ry,Rz):
    roll, pitch = [], []
    for i in range(0,len(Rx)):
        r, p = FSTYAW_1(Rz[i],Rx[i],Ry[i],1)
        roll.append(r)
        pitch.append(p)
        
    return roll,pitch

#execution of the main module

ROLL = []
PITCH = []
YAW = []

ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp4/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp5/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp6/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp7/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp8/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp9/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp10/AQWA004.CSV")))
ROLL.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp11/AQWA004.CSV")))

PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp4/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp5/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp6/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp7/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp8/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp9/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp10/AQWA005.CSV")))
PITCH.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp11/AQWA005.CSV")))

YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp4/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp5/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp6/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp7/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp8/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp9/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp10/AQWA006.CSV")))
YAW.append(readCSV(("E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_submerged/Tp11/AQWA006.CSV")))

#----------------------------------------------------------------------reading the result files 

maxRoll = []
maxPitch = []

roll = []
pitch = []

for i in range(0,len(ROLL)):
    rollT, pitchT = FRAtoLSA(ROLL[i],PITCH[i],YAW[i])
    roll.append(rollT)
    pitch.append(pitchT)
    
for i in range(0,len(roll)):
    maxRoll.append(round(printExtremeValue(roll[i],Globs.FS,Globs.N),2))
    maxPitch.append(round(printExtremeValue(pitch[i],Globs.FS,Globs.N),2))
    
print maxRoll
print maxPitch