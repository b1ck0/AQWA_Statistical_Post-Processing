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
    
def savePDF():
    pp = PdfPages('E:/STORAGE/Python/Statistics/multipage.pdf')
    plt.savefig(pp, format='pdf')
    pp.close()

def MPELine(array,MAX,plotID):
    plt.subplot(plotID)
    x = [0, len(array)]
    y = [MAX, MAX]
    plt.plot(x,y)
    plt.text(len(array)/2.0, 1.02*MAX, 'MPE Value', fontsize=15)

def pltTimeHistoryForce(array,plotID,xLabel,yLabel,title):
    plt.subplot(plotID)
    plt.grid(True)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(title)
    plt.plot(array)
    
def pltTimeHistory(array,plotID,xLabel,yLabel,title):
    plt.subplot(plotID)
    plt.grid(True)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.plot(array)
    
def pltSpectralDensity(array,FS,N,plotID,title,yLabel):
    xf, yf = welch(array,fs=FS,window='hann',nperseg=N,return_onesided=True)
    x,y = fTOw(xf,yf) #angular frequency
    plt.subplot(plotID)
    plt.plot(x, y)
    plt.xlabel('Frequency [rad/s]')
    plt.ylabel(yLabel)
    plt.title(title)
    
def pltExtremeValuesDist(array,FS,N,plotID,yLabel,plotID2,Fig0,Fig1):
    plt.subplot(plotID)
    xf, yf = welch(array,fs=FS,window='hann',nperseg=N,return_onesided=True)
    MAX = Rmax(xf,yf)+np.mean(array)
    p, a = extremeValues(xf, yf)
    plt.plot(p,a+np.mean(array))
    plt.grid(True) 
    plt.ylabel(yLabel)
    plt.xlabel('Probability of Exceedence')
    plt.title('Extreme Values (Rayleigh Distribution)')
    plt.annotate('MPE = ' + str(np.round(MAX,1)), xy=(0.63,MAX), xytext=(0.02 , min(a)+np.mean(array)),arrowprops=dict(facecolor='black', shrink=0.05))
    plt.figure(Fig1)
    MPELine(array,MAX,plotID2)
    plt.figure(Fig0)
        
def pltXY(x,y,plotID,xLabel,yLabel,title):
    plt.subplot(plotID)
    plt.grid(True)
    plt.xlabel(xLabel)
    plt.ylabel(yLabel)
    plt.title(title)
    plt.plot(x,y)
    
def printDispl(x,y,z,Rx,Ry,Rz,filename,figNo,title):
    plt.figure(figNo,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltTimeHistory(x,331,'Time [sec]','X Displacement [m]','X-displacement Time History')
    pltTimeHistory(y,332,'Time [sec]','Y Displacement [m]','Y-displacement Time History')
    pltTimeHistory(z,333,'Time [sec]','Z Displacement [m]','Z-displacement Time History')
    pltTimeHistory(Rx,334,'Time [sec]','Roll Angle [deg]','Roll Angle Time History')
    pltTimeHistory(Ry,335,'Time [sec]','Pitch Angle [deg]','Pitch Angle Time History')
    pltTimeHistory(Rz,336,'Time [sec]','Yaw Angle [deg]','Yaw Angle Time History')
    
    pltXY(x,y,337,'X Displacement [m]', 'Y Displacement [m]', 'Horizontal Position of Ship CoG')
    pltXY(Rz,Rx,338,'Yaw Angle [deg]', 'Roll Angle [deg]', 'Roll Angle vs Yaw Angle')
    pltXY(Rz,Ry,339,'Yaw Angle [deg]', 'Pitch Angle [deg]', 'Pitch Angle vs Yaw Angle')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(figNo)

def printForce(x,y,z,Rx,Ry,Rz,filename,figNo,title):
    plt.figure(figNo,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltTimeHistory(x,331,'Time [sec]','X Load Component [N]','X Load Component Time History')
    pltTimeHistory(y,332,'Time [sec]','Y Load Component [N]','Y Load Component Time History')
    pltTimeHistory(z,333,'Time [sec]','Z Load Component [N]','Z Load Component Time History')
    pltTimeHistory(Rx,334,'Time [sec]','Mx Load Component [N.m]','Mx Load Component Time History')
    pltTimeHistory(Ry,335,'Time [sec]','My Load Component [N.m]','My Load Component Time History')
    pltTimeHistory(Rz,336,'Time [sec]','Mz Load Component [N.m]','Mz Load Component Time History')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(figNo)
    
def printTugForce(Fx,Fy,filename,figNo,title):
    plt.figure(figNo,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltTimeHistory(Fx,311,'Time [sec]','X Force [N]','X Force Time History')
    pltTimeHistory(Fy,312,'Time [sec]','Y Force [N]','Y Force Time History')
   
    F = []
    
    for i in range(0,len(Fx)):
        F.append(sqrt(Fx[i]**2+Fy[i]**2))

    pltTimeHistory(F,313,'Time [sec]','Total Force [N]','Total Force Time History')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(figNo)

def printSlings(gTen1,gTen2,filename,figNo,title):
    plt.figure(figNo,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')

    pltTimeHistoryForce(gTen1,321,'Time [sec]','FWD Hook Load [N]','FWD Hook Load Time Hostory')
    pltTimeHistoryForce(gTen2,322,'Time [sec]','AFT Hook Load [N]','AFT Hook Load Time Hostory')
    
    pltSpectralDensity(gTen1,Globs.FS,Globs.N,323,'Power Spectral Density FWD Hook Load','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(gTen2,Globs.FS,Globs.N,324,'Power Spectral Density AFT Hook Load','Power Spectral Density [N$^2$/s]')
    
    pltExtremeValuesDist(gTen1,Globs.FS,Globs.N,325,'FWD Hook Load [N]',321,figNo,figNo)
    pltExtremeValuesDist(gTen2,Globs.FS,Globs.N,326,'FWD Hook Load [N]',322,figNo,figNo)
        
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(filename,bbox_inches='tight')
    plt.close(figNo)
    
def printSlingsCPlate(gTen1,gTen2,gTen3,gTen4,directory,figNo1,figNo2,figNo3,title):
    plt.figure(figNo1,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')

    pltTimeHistoryForce(gTen1,221,'Time [sec]','FWD Hook Load [N]','FWD/STB Load Time Hostory')
    pltTimeHistoryForce(gTen2,222,'Time [sec]','AFT Hook Load [N]','FWD/PS Hook Load Time Hostory')
    pltTimeHistoryForce(gTen3,223,'Time [sec]','AFT Hook Load [N]','AFT/STB Hook Load Time Hostory')
    pltTimeHistoryForce(gTen4,224,'Time [sec]','AFT Hook Load [N]','AFT/PS Hook Load Time Hostory')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/09_cPlate_slings_TH.png',bbox_inches='tight')
    
    plt.figure(figNo2,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltSpectralDensity(gTen1,Globs.FS,Globs.N,221,'Power Spectral Density FWD/STB Hook Load','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(gTen2,Globs.FS,Globs.N,222,'Power Spectral Density FWD/PS Hook Load','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(gTen3,Globs.FS,Globs.N,223,'Power Spectral Density AFT/STB Hook Load','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(gTen4,Globs.FS,Globs.N,224,'Power Spectral Density AFT/PS Hook Load','Power Spectral Density [N$^2$/s]')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/10_cPlate_slings_FD.png',bbox_inches='tight')

    plt.figure(figNo3,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltExtremeValuesDist(gTen1,Globs.FS,Globs.N,221,'FWD/STB Hook Load [N]',221,figNo3,figNo1)
    pltExtremeValuesDist(gTen2,Globs.FS,Globs.N,222,'FWD/PS Hook Load [N]',222,figNo3,figNo1)
    pltExtremeValuesDist(gTen3,Globs.FS,Globs.N,223,'AFT/STB Hook Load [N]',223,figNo3,figNo1)
    pltExtremeValuesDist(gTen4,Globs.FS,Globs.N,224,'AFT/PS Hook Load [N]',224,figNo3,figNo1)
        
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/11_cPlate_slings_STAT.png',bbox_inches='tight')
    
    plt.close(figNo1)
    plt.close(figNo2)
    plt.close(figNo3)
     
def printMoorings(ten0,ten1,uplift0,uplift1,llength0,llength1,directory,figNo1,figNo2,figNo3,title):
    plt.figure(figNo1,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltTimeHistoryForce(ten0, 321, 'Time [sec]', 'Tension [N]','STB Mooring Line Time History')
    pltTimeHistoryForce(ten1, 322, 'Time [sec]', 'Tension [N]','PS Mooring Line Time History')
    pltTimeHistoryForce(uplift0, 323, 'Time [sec]', 'Anchor Uplift [N]','STB Mooring Line Time History')
    pltTimeHistoryForce(uplift1, 324, 'Time [sec]', 'Anchor Uplift [N]','PS Mooring Line Time History')
    
    pltTimeHistory(llength0, 325, 'Time [sec]', 'Laid Length [m]','STB Mooring Line Time History')
    pltTimeHistory(llength1, 326, 'Time [sec]', 'Laid Length [m]','PS Mooring Line Time History')   
        
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/05_ship_mooring_TH.png',bbox_inches='tight')

    plt.figure(figNo2,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltSpectralDensity(ten0,Globs.FS,Globs.N,321,'Power Spectral Density STB Mooring Line Tension','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(ten1,Globs.FS,Globs.N,322,'Power Spectral Density PS Mooring Line Tension','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(uplift0,Globs.FS,Globs.N,323,'Power Spectral Density STB Anchor Uplift','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(uplift1,Globs.FS,Globs.N,324,'Power Spectral Density PS Anchor Uplift','Power Spectral Density [N$^2$/s]')
    pltSpectralDensity(llength0,Globs.FS,Globs.N,325,'Power Spectral Density PS Laid Length','Power Spectral Density [m$^2$/s]')
    pltSpectralDensity(llength1,Globs.FS,Globs.N,326,'Power Spectral Density STB Laid Length','Power Spectral Density [m$^2$/s]')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/06_ship_mooring_FD.png',bbox_inches='tight')
    
    plt.figure(figNo3,figsize=(25,12))
    plt.suptitle(title,fontsize=14,fontweight='bold')
    
    pltExtremeValuesDist(ten0,Globs.FS,Globs.N,321,'Tension STB Mooring Line [N]',321,figNo3,figNo1)
    pltExtremeValuesDist(ten1,Globs.FS,Globs.N,322,'Tension PS Mooring Line [N]',322,figNo3,figNo1)
    pltExtremeValuesDist(uplift0,Globs.FS,Globs.N,323,'Anchor Uplift STB Mooring Line [N]',323,figNo3,figNo1)
    pltExtremeValuesDist(uplift1,Globs.FS,Globs.N,324,'Anchor Uplift PS Mooring Line [N]',324,figNo3,figNo1)
    pltExtremeValuesDist([Globs.L - i for i  in llength0],Globs.FS,Globs.N,325,'Lifted Length STB Mooring Line [m]',325,figNo3,figNo1)
    pltExtremeValuesDist([Globs.L - i for i  in llength1],Globs.FS,Globs.N,326,'Lifted Length PS Mooring Line [m]',326,figNo3,figNo1)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(directory + '/07_ship_mooring_STAT.png',bbox_inches='tight')
    
    plt.close(figNo1)
    plt.close(figNo2)
    plt.close(figNo3)
        
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
    
def extremeValues(xf, yf):
    a = []
    prob = np.arange(0.01,1.0,0.01) #probabilities for not exceedence
    t = Globs.t
    yf2 = yf * xf**2 #ordinates of the second order spectrum
    M0 = np.trapz(yf,xf) #integrating the area under the yf curve
    M2 = np.trapz(yf2,xf) #integrating the area under the yf2 curve
    std = np.sqrt(M0) #obtaining the standard deviation of the input signal = square root of the variance
    Vx = np.sqrt(M2/M0) #average zero-upcrossing frequency of the input signal
    for p in prob:
        a.append(std*np.sqrt(2.0*np.log((Vx*t)/(np.log(1.0/p))))) 
        
    return prob[::-1],a

def fTOw(xf,yf):
    x = xf * 2.0*np.pi
    y = yf / (2.0*np.pi)
    
    return x,y
       
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
    for i in range(0,len(Rx[0])):
        r, p = FSTYAW_1(Rz[0][i],Rx[0][i],Ry[0][i],1)
        roll.append(r)
        pitch.append(p)
        
    return roll,pitch

class Ship:
    def __init__(self):
        self.x = [] #container for X position
        self.y = [] #container for Y position
        self.z = [] #container for Z position
        self.Rx = [] #container for Rx position
        self.Ry = [] #container for Ry position
        self.Rz = [] #container for Rz position
        self.mTenSTB = [] #container for Mooring Tension 1 (STB Line)
        self.mTenPS = [] #container for Mooring Tension 2 (PS Line)
        self.mUpliftSTB = [] #container for Mooring Laid Length 1 (STB Line)
        self.mUpliftPS = [] #container for Mooring Laid Length 2 (PS Line)
        self.mLLengthSTB = [] #container for Mooring Laid Length 1 (STB Line)
        self.mLLengthPS = [] #container for Mooring Laid Length 2 (PS Line)
        self.gTenFWD = [] #container for Lifting Grommet Tension (FWD Line)
        self.gTenAFT = [] #container for Lifting Grommet Tension (AFT Line)
        self.tTen = [] #container for Tugline Tension
        self.p1X = [] #container for X-coordinate of Crane Point 1 in FRA
        self.p1Y = [] #container for Y-coordinate of Crane Point 1 in FRA
        self.p2X = [] #container for X-coordinate of Crane Point 2 in FRA
        self.p2Y = [] #container for Y-coordinate of Crane Point 2 in FRA
        self.roll = []
        self.pitch = []
           
class Tug:
    def __init__(self):
        self.x = [] #container for X position
        self.y = [] #container for Y position
        self.z = [] #container for Z position
        self.Rx = [] #container for Rx position
        self.Ry = [] #container for Ry position
        self.Rz = [] #container for Rz position
        self.Fx = [] #container for X-component of the External Load
        self.Fy = [] #container for Y-component of the External Load
        
class CoverPlate:
    def __init__(self):
        self.x = [] #container for X position
        self.y = [] #container for Y position
        self.z = [] #container for Z position
        self.Rx = [] #container for Rx position
        self.Ry = [] #container for Ry position
        self.Rz = [] #container for Rz position
        self.MorX = []
        self.MorY = []
        self.MorZ = []
        self.MorRx = []
        self.MorRy = []
        self.MorRz = []
        self.SlamX = []
        self.SlamY = []
        self.SlamZ = []
        self.SlamRx = []
        self.SlamRy = []
        self.SlamRz = []
        self.gTen1 = []
        self.gTen2 = []
        self.gTen3 = []
        self.gTen4 = []
        self.p1X = []
        self.p1Y = []
        self.p2X = []
        self.p2Y = []


#execution of the main module

ship = Ship()
tug = Tug()
cPlate = CoverPlate()

directory = "E:/STORAGE/ANSYS PROJECTS/H5_0117_Husky_R1/results_splashZone/Tp11"

#----------------------------------------------------------------------reading the result files 

#----------------------------------------reading the files for the Ship

ship.x.append(readCSV(directory + "/AQWA001.CSV")) # X-coordinate of CoG
ship.y.append(readCSV(directory + "/AQWA002.CSV")) # Y-coordinate of CoG
ship.z.append(readCSV(directory + "/AQWA003.CSV")) # Z-coordinate of CoG
ship.Rx.append(readCSV(directory + "/AQWA004.CSV")) # Rx-angle in FRA
ship.Ry.append(readCSV(directory + "/AQWA005.CSV")) # Ry-angle in FRA
ship.Rz.append(readCSV(directory + "/AQWA006.CSV")) # Rz-angle in FRA
ship.mTenSTB.append(readCSV(directory + "/AQWA007.CSV")) # Mooring Tension 1 (STB Line)
ship.mUpliftSTB.append(readCSV(directory + "/AQWA008.CSV")) # Mooring Uplift 1 (STB Line)
ship.mLLengthSTB.append(readCSV(directory + "/AQWA009.CSV")) # Mooring Laid Length 1 (STB Line)
ship.mTenPS.append(readCSV(directory + "/AQWA010.CSV")) # Mooring Tension 2 (PS Line)
ship.mUpliftPS.append(readCSV(directory + "/AQWA011.CSV")) # Mooring Uplift 2 (PS Line)
ship.mLLengthPS.append(readCSV(directory + "/AQWA012.CSV")) # Mooring Laid Length 2 (PS Line)
ship.gTenFWD.append(readCSV(directory + "/AQWA013.CSV")) # Lifting Grommet Tension (FWD Line)
ship.gTenAFT.append(readCSV(directory + "/AQWA014.CSV")) # Lifting Grommet Tension (AFT Line)
ship.tTen.append(readCSV(directory + "/AQWA015.CSV")) # Tugline Tension
# X-coordinate of Crane Point 1 in FRA
# Y-coordinate of Crane Point 1 in FRA
# X-coordinate of Crane Point 2 in FRA
# Y-coordinate of Crane Point 2 in FRA


#----------------------------------------reading the files for the Cover Plate

cPlate.x.append(readCSV(directory + "/AQWA016.CSV")) # X-coordinate of CoG
cPlate.y.append(readCSV(directory + "/AQWA017.CSV")) # Y-coordinate of CoG
cPlate.z.append(readCSV(directory + "/AQWA018.CSV")) # Z-coordinate of CoG
cPlate.Rx.append(readCSV(directory + "/AQWA019.CSV")) # Rx-angle in FRA
cPlate.Ry.append(readCSV(directory + "/AQWA020.CSV")) # Ry-angle in FRA
cPlate.Rz.append(readCSV(directory + "/AQWA021.CSV")) # Rz-angle in FRA
cPlate.MorX.append(readCSV(directory + "/AQWA022.CSV")) # X-component of Morison Drag 
cPlate.MorY.append(readCSV(directory + "/AQWA023.CSV")) # Y-component of Morison Drag 
cPlate.MorZ.append(readCSV(directory + "/AQWA024.CSV")) # Z-component of Morison Drag 
cPlate.MorRx.append(readCSV(directory + "/AQWA025.CSV")) # Rx-component of Morison Drag 
cPlate.MorRy.append(readCSV(directory + "/AQWA026.CSV")) # Ry-component of Morison Drag
cPlate.MorRz.append(readCSV(directory + "/AQWA027.CSV")) # Rz-component of Morsion Drag
cPlate.SlamX.append(readCSV(directory + "/AQWA028.CSV")) # X-component of Slam 
cPlate.SlamY.append(readCSV(directory + "/AQWA029.CSV")) # Y-component of Slam 
cPlate.SlamZ.append(readCSV(directory + "/AQWA030.CSV")) # Z-component of Slam 
cPlate.SlamRx.append(readCSV(directory + "/AQWA031.CSV")) # Rx-component of Slam
cPlate.SlamRy.append(readCSV(directory + "/AQWA032.CSV")) # Ry-component of Slam
cPlate.SlamRz.append(readCSV(directory + "/AQWA033.CSV")) # Rz component of Slam
cPlate.gTen1.append(readCSV(directory + "/AQWA034.CSV")) # Lifting Grommet Tension 1
cPlate.gTen2.append(readCSV(directory + "/AQWA035.CSV")) # Lifting Grommet Tension 2
cPlate.gTen3.append(readCSV(directory + "/AQWA036.CSV")) # Lifting Grommet Tension 3
cPlate.gTen4.append(readCSV(directory + "/AQWA037.CSV")) # Lifting Grommet Tension 4
# X-coordinate of FWD Point in FRA
# Y-coordinate of FWD Point in FRA
# X-coordinate of AFT Point in FRA
# Y-coordinate of AFT Point in FRA

#----------------------------------------reading the files for the TUG

tug.x.append(readCSV(directory + "/AQWA038.CSV")) # X-coordinate of CoG
tug.y.append(readCSV(directory + "/AQWA039.CSV")) # Y-coordinate of CoG
tug.z.append(readCSV(directory + "/AQWA040.CSV")) # Z-coordinate of CoG
tug.Rx.append(readCSV(directory + "/AQWA041.CSV")) # Rx-angle in FRA
tug.Ry.append(readCSV(directory + "/AQWA042.CSV")) # Ry-angle in FRA
tug.Rz.append(readCSV(directory + "/AQWA043.CSV")) # Rz-angle in FRA
tug.Fx.append(readCSV(directory + "/AQWA044.CSV")) # X-component of the External Load
tug.Fy.append(readCSV(directory + "/AQWA045.CSV")) # Y-component of the External Load

#generating the graphical data

filename = directory + "/01_ship_displacement.png"
roll,pitch = FRAtoLSA(ship.Rx,ship.Ry,ship.Rz) #conversion from FRA to LSA
printDispl(ship.x[0],ship.y[0],ship.z[0],roll,pitch,ship.Rz[0],filename,1,"Vessel Displacement Time History")

filename = directory + "/02_cPlate_displacement.png"
roll,pitch = FRAtoLSA(cPlate.Rx,cPlate.Ry,cPlate.Rz) #conversion from FRA to LSA
printDispl(cPlate.x[0],cPlate.y[0],cPlate.z[0],roll,pitch,cPlate.Rz[0],filename,2,"Cover Plate Displacement Time History")

filename = directory + "/03_tug_displacement.png"
roll,pitch = FRAtoLSA(tug.Rx,tug.Ry,tug.Rz) #conversion from FRA to LSA
printDispl(tug.x[0],tug.y[0],tug.z[0],roll,pitch,tug.Rz[0],filename,3,"Tugboat Displacement Time History")

filename = directory + "/04_tug_force.png"
printTugForce(tug.Fx[0],tug.Fy[0],filename,4,"External Forces on the Tugboat Time History")

printMoorings(ship.mTenSTB[0],ship.mTenPS[0],ship.mUpliftSTB[0],ship.mUpliftPS[0],ship.mLLengthSTB[0],ship.mLLengthPS[0],directory,5,6,7,"Vessel Mooring Forces")

filename = directory + "/08_ship_slings.png"
printSlings(ship.gTenFWD[0],ship.gTenAFT[0],filename,8,"Steel Grommets Forces")

printSlingsCPlate(cPlate.gTen1[0],cPlate.gTen2[0],cPlate.gTen4[0],cPlate.gTen4[0],directory,9,10,11, "Polyester Grommets Forces")

filename = directory + "/12_cPlate_Morison.png"
printForce(cPlate.MorX[0],cPlate.MorY[0],cPlate.MorZ[0],cPlate.MorRx[0],cPlate.MorRy[0],cPlate.MorRz[0],filename,12,"Morison Forces on the Cover Plate Time History")

filename = directory + "/13_cPlate_Slam.png"
printForce(cPlate.SlamX[0],cPlate.SlamY[0],cPlate.SlamZ[0],cPlate.SlamRx[0],cPlate.SlamRy[0],cPlate.SlamRz[0],filename,13, "Slamming Forces on the Cover Plate Time History")