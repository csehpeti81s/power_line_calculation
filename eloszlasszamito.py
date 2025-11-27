#!/usr/bin/env python3 
import math, cmath, numpy
import pandas as texout

DEBUG=False
ONE_LINE_CALCULATION = False

def replusZ(Z1, Z2):
    return (Z1*Z2)/(Z1+Z2)

def GMRfromA(A):
    return 0.78*cmath.sqrt(A/1000/1000/cmath.pi)

def selfReactance(D, GMR):
    return 0.145*cmath.log10(D/GMR)

def mutReactance(D_e, D_fg):
    return 0.145*cmath.log10(D_e/D_fg)

def currentRatio(Z1, Z2):
    Ir1 = Z2 / (Z1+Z2)
    Ir2 = Z1 / (Z1+Z2)
    return (Ir1, Ir2)

def genuineRoot(Z1, Z2):
    if (Z1.real >= 0) & (Z2.real >= 0):
        raise "Mindkét gyök valós!"
    if Z1.real >= 0:
        return Z1
    if Z2.real >=0:
        return Z2
    raise "Nincs valós ellenállást tartalmazó gyök"

def quadraticFormula(a, b, c):
    D= (b*b-4*a*c)
    x1=(-b + cmath.sqrt(D))/(2*a)
    x2=(-b - cmath.sqrt(D))/(2*a)
    return (x1, x2)

class cPowerLine():
    def __init__(self, towerCount):
        self.U = 76.21
        self.Izsubst = complex(0, -31.5)
        self.Izref = self.Izsubst
        self.Zmh = self.U /self.Izsubst
        self.Uphwline = None
        self.Ustwline = None
        self.Uz = None
        self.Zz = None
        self.Iz = self.Izref
        self.lineLength = 0
        self.towerList = None
        self.Igmax = 0
        self.Igmaxtwrid = 0
        self.makePowerLine(towerCount)
        self.updateTowerInternalCalculatedValues()
        for i in range(5):
            self.updatePowerLine()
            if DEBUG:
                dbgtxt = (f'{"*" * 30}{i:^5d}{"*" * 30}')
                print(dbgtxt)
                self.printIntermediateResultsOfTowers()
    def updatePowerLine(self):
        self.updateShortCircuitValuesToIz()
        self.calculateUeZe()
        self.calculateIgrw()
        self.calculateLoopVoltages()
        self.calculateShortCircuitCurrent()
        self.calculateLineValues()
        if DEBUG:
            dbgmsg = (f'{"N":>5} {"Uz":>15} {"Zz":>15} {"Izprev":>15} {"Iz":>15} \n'
                f'{len(self.towerList):5d} '
                f'{abs(self.Uz):9.2f}<{cmath.phase(self.Uz):5.2f} {self.Zz:15.2f} '
                f'{self.towerList[0].Izref:15.2f} {self.Iz:15.2f}')
            print(dbgmsg)
    def makePowerLine(self, towerCount):
        self.towerList = []              #Towers and spans are counted from the faulty ones. 
        for i in range(towerCount):
            tower = cTower()
            tower.Izref = self.Izref
            self.towerList.append(tower) 
        self.towerList[-1].Dtowers = 0.10  #Fist tower is closer to substation than others. 
        #self.towerList[0].Dgrr = self.towerList[0].Dgrr / 2 #Depth is half at tower 0. 
    def updateTowerInternalCalculatedValues(self):
        for tower in self.towerList:
            tower.calcInternalValues()
    def calculateUeZe(self):
        nextTower = None
        for tower in reversed(self.towerList):
            tower.calculateZeUe(nextTower)
            nextTower = tower
    def calculateIgrw(self):
        previousTower = None
        for tower in self.towerList:
            tower.calculateStaticWireCurrent(previousTower)
            previousTower = tower
    def calculateLoopVoltages(self):
        self.Uz = self.towerList[-1].Istwsum * self.towerList[-1].Rsubstgr
        for tower in self.towerList:
            tower.calculateSpanVoltages()
            self.Uz =self.Uz + tower.Uspan
    def calculateShortCircuitCurrent(self):
        self.Zz = self.Uz / self.Iz
        self.Iz = self.U / (self.Zz + self.Zmh )
    def calculateLineValues(self):
        self.lineLength = 0
        self.Igmax = 0
        for idx, tower in enumerate(self.towerList):
            self.lineLength = self.lineLength + int(tower.Dtowers *1000)
            if (abs(self.Igmax) < abs(tower.Istwsum)):
                self.Igmax = abs(tower.Istwsum)
                self.Igmaxtwrid = idx
    def updateShortCircuitValuesToIz(self):
        for tower in self.towerList:
            tower.Izref = self.Iz
    def printIntermediateResultsOfTowers(self):
        headertxt = (f'{"#":>3}{"Táv":>4}'
            f'{"Zstw":>13}{"Zphw":>13}'
            f'{"Zmut":>13}{"Ustphw":>13}'
            f'{"Ze":>13}{"Ue":>13}'
            f'{"Istwind":>13}{"Istwcon":>13}{"Ist":>13}'
            f'{"Uphstw":>16}'
            f'{"Uspan":>16}{"Izreftwr":>16}')
        msgtxt = headertxt + "\n"
        for idx, tower in enumerate(self.towerList):
            towertxt = (f'{idx:3d}{tower.Dtowers*1000:4.0f}'
                    f'{tower.Zstw:13.2f}{tower.Zphw:13.2f}' 
                    f'{tower.Zphstmut:13.2f}{(tower.Ustphw):13.2f}'
                    f'{(tower.Ze):13.2f}{(tower.Ue):13.2f}'
                    f'{tower.Istwind:13.2f}{tower.Istwcond:13.2f}'
                    f'{(tower.Istwind+tower.Istwcond):13.2f}'
                    f'{(tower.Uphstw):16.2f}'
                    f'{tower.Uspan:16.2f}{tower.Izref:16.2f}' )
            msgtxt = msgtxt + towertxt + "\n"
        with open("towermsg.txt", "w") as f:
            f.write(msgtxt)
    def printIntermediateResultsOfLine(self):
        msgtxt = (f'{"Twr c":>6}{"Izref":>13}'
            f'{"Uz":>22}{"Zz":>22}{"Iz":>22}\n')
        msgtxt = msgtxt + (
            f'{len(self.towerList):6}{self.Izref:13.2f}'
            f'{self.Uz:16.3f}({abs(self.Uz):4.2f}){self.Zz:16.2f}({abs(self.Zz):4.2f})'
            f'{self.Iz:16.2f}({abs(self.Iz):4.2f})\n') 
        with open("linemsg.txt", "w") as f:
            f.write(msgtxt)

class cTower():
    def __init__(self):
        self.Izref = None         #reference value of earth fault current in kA 
        self.Dtowers = 0.3    #lengt in kilometers between this and next pole
        self.Rtowergr = 5        #earth resistance of the current pole in ohm
        self.Rsubstgr = 0.05    #earth resistance of the substation in ohm
        self.Astw = 95 + 50     #cross section of the static wire in mm2
        self.Aphw = 250 + 40    #cross section of the phase wire in mm2
        self.Rphw = 0.1154      #resistivity of the phase wire in ohm/km
        self.Rstw = 0.2992      #resistivity of the static wire in ohm/km
        self.Rgrr = 0           #resistivity of earth return in ohm/km
        self.Dphstw = 10        #distance of ground and static wire in meter
        self.Dgrr = 1000        #depth of earth return current in meter
        self.Zstw = None        #self impedance of static wire ground loop ohm/span
        self.Zphw = None        #sefl ipmedance of phase wire ground loop ohm/span
        self.Zphstmut = None    #mutual impedance of phase and static to ground loops ohm/span 
        self.Ustphw = None
        self.Uphstw = None
        self.Ue = None 
        self.Ze = None
        self.Istwind = None
        self.Istwcond = None
        self.Istwsum = None
        self.Ustwspan = None        #Sum voltage on static wire
        self.Uphwspan = None
        self.Uspan = None
    def calcInternalValues(self):
        self.Zstw =  self.Dtowers * complex(
            self.Rstw, selfReactance(self.Dgrr, GMRfromA(self.Astw)))  
        self.Zphw = self.Dtowers * complex(
            self.Rphw, selfReactance(self.Dgrr, GMRfromA(self.Aphw))) 
        self.Zphstmut = self.Dtowers * complex(
            self.Rgrr, mutReactance(self.Dgrr, self.Dphstw))
    def calculateZeUe(self, next):
        self.Ustphw = self.Zphstmut * self.Izref
        if next is None:    #This case I'am next to the substation.
            self.Ze = self.Zstw + self.Rsubstgr  
            self.Ue = self.Ustphw
            return
        self.Ze = self.Zstw + replusZ(self.Rtowergr, next.Ze)
        self.Ue = self.Ustphw + next.Ue * self.Rtowergr / (self.Rtowergr + next.Ze)
    def calculateStaticWireCurrent(self, previous):
        currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        if DEBUG:
            dbgmsg = f'Current ratio in static wire: {currRatio:5.2f}'
            print(dbgmsg)
        if previous is None:    #Means I'am at the faulty tower. 
            self.Istwcond = self.Izref * currRatio
            self.Istwind = self.Ue / (self.Rtowergr + self.Ze)
            self.Istwsum = (self.Istwcond + self.Istwind)
            self.Uphstw = self.Zphstmut * self.Istwsum
            return
        self.Istwcond = previous.Istwcond * currRatio
        #currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        self.Istwind = currRatio * previous.Istwind + self.Ue / (self.Rtowergr + self.Ze)
        self.Istwsum = (self.Istwcond + self.Istwind)
        self.Uphstw = self.Zphstmut * (self.Istwsum)
    def calculateSpanVoltages(self):
        #positive values are in the Iz direction
        self.Ustwspan = self.Istwsum * self.Zstw - self.Izref * self.Zphstmut
        self.Uphwspan = self.Izref * self.Zphw - self.Istwsum * self.Zphstmut
        self.Uspan = self.Uphwspan + self.Ustwspan
        if DEBUG:
            dbgtxt = (f'Ustwspan = Istw * Zstw - Izref * Zphstmut \n'
                f'{self.Ustwspan:.2f}={self.Istwsum:.2f}*{self.Zstw:.2f}'
                f'-{self.Izref:.2f}*{self.Zphstmut:.2f}\n'
                f'Uphwspan = Izref * Zphw - Istw * Zphstmut \n'
                f'{self.Uphwspan:.2f}={self.Izref:.2f}*{self.Zphw:.2f}'
                f'-{self.Istwsum:.2f}*{self.Zphstmut:.2f}\n'
                f'Uspan={self.Uspan:.2} |Uspan|={abs(self.Uspan)}\n'
                f'{abs(self.Uphwspan/self.Uspan):5.2f}'
                f'{abs(self.Ustwspan/self.Uspan):5.2f}'
                f'{abs(self.Ustwspan/self.Izref):5.2f}'
                f'{abs(self.Uphwspan/self.Izref):5.2f}'
                f'{abs((self.Uphwspan+self.Ustwspan)/self.Izref):5.2f}'
                )
            print(dbgtxt)

def exportTowerResultsToLatex(towerList):
    towerid = []
    towerdistance = []
    Zphw = []
    Zstw = []
    Zstphw = []
    Istwind = []
    Istwcond = []
    Istwsum = []
    for idx, tower in enumerate(towerList):
        towerid.append(idx)
        towerdistance.append(int(tower.Dtowers*1000))
        Zphw.append(abs(tower.Zphw))
        Zstw.append(abs(tower.Zstw))
        Zstphw.append(abs(tower.Zphstmut))
        Istwind.append(abs(tower.Istwind))
        Istwcond.append(abs(tower.Istwcond))
        Istwsum.append(abs(tower.Istwsum))
    data = {
            "Oszlop": tuple(towerid),
            "Oszlopköz": tuple(towerdistance),
            "Zphw": tuple(Zphw),
            "Zstw": tuple(Zstw),
            "Zstphw": tuple(Zstphw),
            "Ig kond.": tuple(Istwcond),
            "Ig indukt.": tuple(Istwind),
            "Ig": tuple(Istwsum),
        }
    filename = (f'./out/tavvezetek{len(towerList):03d}py.tex')
    writeToLatex(data, filename)

def writeToLatex(data, filename):
    formatteddata = texout.DataFrame(data)
    latextxt = formatteddata.to_latex(
        index = False,
        position = "h",
        float_format="%.2f",
        #caption="Távvezeték zárlati áramok nagysága és eloszlása",
        #label="tab:tavvezetek_eredmenyek"
        ).replace(".", ",")
    print(latextxt)
    with open(filename, "w") as f:
        f.write(latextxt)
    
def calculateSinglePowerLine(towerCount):
    powerLine = cPowerLine(towerCount)
    powerLine.printIntermediateResultsOfTowers()
    powerLine.printIntermediateResultsOfLine()
    exportTowerResultsToLatex(powerLine.towerList)

def calculateMultiplePowerLines():
    towerSet = tuple(range(25)) + (49, 99, 199)    #.append(100).append(200) 
    towerCounts = []
    lineLengths = []
    Uz = []
    Igmax = []
    Igmaxtwrid = []
    Iz = []
    Izrel = []
    for towerCount in towerSet:
        powerLine = cPowerLine(towerCount+1)
        towerCounts.append(len(powerLine.towerList))
        lineLengths.append(powerLine.lineLength)
        Igmax.append(powerLine.Igmax)
        Igmaxtwrid.append(powerLine.Igmaxtwrid)
        Uz.append(abs(powerLine.Uz))
        Iz.append(abs(powerLine.Iz))
        Izrel.append(abs(powerLine.Iz) / abs(powerLine.Izsubst))
    data = {
        "Oszlopszám" : tuple(towerCounts),
        "Vonalhossz": tuple(lineLengths),
        "Igmax": tuple(Igmax),
        "Igmax oszl": tuple(Igmaxtwrid),
        "Uz": tuple(Uz),
        "Iz": tuple(Iz),
        "Izrel": tuple(Izrel),
    }
    filename = f'./out/sorozat.tex'
    writeToLatex(data, filename)
    
def main():
    if ONE_LINE_CALCULATION:
        towerset = [3,5,10,50]
        if DEBUG:
            towerset = []
            towerset.append(2)
        for towerCount in (towerset):
            calculateSinglePowerLine(towerCount)
        return    
    calculateMultiplePowerLines()
    

if __name__=="__main__":
    main()
