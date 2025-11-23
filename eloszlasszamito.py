#!/usr/bin/env python3 
import math, cmath, numpy
import pandas as texout

DEBUG=True
#DEBUG=False

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
        self.Iz = None
        self.towerList = None
        self.makePowerLine(towerCount)
        self.updateTowerInternalCalculatedValues()
        self.calculateUeZe()
        self.calculateIgrw()
        self.calculateLoopVoltages()
        self.calculateShortCircuitCurrent()
        self.updateShortCircuitValuesToIz()
        self.updateTowerInternalCalculatedValues()
        self.calculateUeZe()
        self.calculateIgrw()
        self.calculateLoopVoltages()
        self.calculateShortCircuitCurrent()
    def makePowerLine(self, towerCount):
        self.towerList = []              #Towers and spans are counted from the faulty ones. 
        for i in range(towerCount):
            tower = cTower()
            tower.Izref = self.Izref
            self.towerList.append(tower) 
        self.towerList[-1].Dtowers = 0.1  #Fist tower is closer to substation than others. 
        self.towerList[0].Dgrr = self.towerList[0].Dgrr / 2 #Depth is half in the faulty tower span. 
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
            tower.calculateShWCurrent(previousTower)
            previousTower = tower
    def calculateLoopVoltages(self):
        for tower in self.towerList:
            tower.calculateSpanVoltages()
        self.Uphwline = complex(0,0)
        self.Ustwline = complex(0,0)
        for tower in self.towerList:
            self.Uphwline = self.Uphwline + tower.Uphwspan
            self.Ustwline = self.Ustwline + tower.Ustwspan
        self.Uz = self.Uphwline + self.Ustwline
    def calculateShortCircuitCurrent(self):
        self.Zz = self.Uz / self.Izref
        self.Iz = self.U / (self.Zz + self.Zmh )
    def updateShortCircuitValuesToIz(self):
        for tower in self.towerList:
            tower.Izref = self.Iz
    def printIntermediateResults(self):
        headertxt = (f'{"#":>3}{"Táv":>4}'
            f'{"Zstw":>13}{"Zphw":>13}'
            f'{"Zmut":>13}{"Ustphw/k":>13}'
            f'{"Ze":>13}{"Ue/k":>13}'
            f'{"Istwmut/k":>13}{"Istwcon/k":>13}{"Ist/k":>13}'
            f'{"Uphstw/k":>16}{"Uphwspan/k":>16}{"Ustwspan/k":>16}')
        msgtxt = headertxt + "\n"
        for idx, tower in enumerate(self.towerList):
            towertxt = (f'{idx:3d}{tower.Dtowers*1000:4.0f}'
                    f'{tower.Zstw:13.2f}{tower.Zphw:13.2f}' 
                    f'{tower.Zphstmut:13.2f}{(tower.Ustphw):13.2f}'
                    f'{(tower.Ze):13.2f}{(tower.Ue):13.2f}'
                    f'{tower.Istwmut:13.2f}{tower.Istwcond:13.2f}{(tower.Istwmut+tower.Istwcond):13.2f}'
                    f'{(tower.Uphstw):16.2f}'
                    f'{tower.Uphwspan:16.2f}{tower.Ustwspan:16.2f}' )
            msgtxt = msgtxt + towertxt + "\n"
        with open("towermsg.txt", "w") as f:
            f.write(msgtxt)
    def printIntermediateResultsOfLine(self):
        msgtxt = (f'{"Twr c":>6}{"Iz/k":>13}{"Ustwline":>13}{"Uphwline":>16}'
            f'{"Uz":>22}{"Zz":>22}{"Iz":>22}\n')
        msgtxt = msgtxt + (
            f'{len(self.towerList):6}{self.Izref:13.2f}{self.Ustwline:13.3f}'
            f'{self.Uphwline:16.3f}'
            f'{self.Uz:16.3f}({abs(self.Uz):4.2f}){self.Zz:16.2f}({abs(self.Zz):4.2f})'
            f'{self.Iz:16.2f}({abs(self.Iz):4.2f})\n') 
        with open("linemsg.txt", "w") as f:
            f.write(msgtxt)

class cTower():
    def __init__(self):
        self.Izref = None         #reference value of earth fault current in kA 
        self.Dtowers = 0.300     #lengt in kilometers between this and previous poles
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
        self.Istwmut = None
        self.Istwcond = None
        self.Istwsum = None
        self.Ustwspan = None        #Sum voltage on static wire
        self.Uphwspan = None
    def calcInternalValues(self):
        self.Zstw =  self.Dtowers * complex(
            self.Rstw, selfReactance(self.Dgrr, GMRfromA(self.Astw)))  
        self.Zphw = self.Dtowers * complex(
            self.Rphw, selfReactance(self.Dgrr, GMRfromA(self.Aphw))) 
        self.Zphstmut = self.Dtowers * complex(
            self.Rgrr, mutReactance(self.Dgrr, self.Dphstw))
        self.Ustphw = self.Zphstmut * self.Izref
    def calculateZeUe(self, next):
        if next is None:
            self.Ze = self.Zstw + self.Rsubstgr + self.Rgrr 
            self.Ue = self.Ustphw
            return
        self.Ze = self.Zstw + self.Rsubstgr + replusZ(self.Rtowergr, next.Ze)
        self.Ue = self.Ustphw + next.Ue * self.Rtowergr / (self.Rtowergr + next.Ze)
    def calculateShWCurrent(self, previous):
        currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        if previous is None:
            self.Istwcond = self.Izref * currRatio
            self.Istwmut = self.Ue / self.Rtowergr
            self.Istwsum = (self.Istwcond + self.Istwmut)
            self.Uphstw = self.Zphstmut * (self.Istwcond + self.Istwmut)
            return
        self.Istwcond = previous.Istwcond * currRatio
        currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        self.Istwmut = currRatio * previous.Istwmut + self.Ue / (self.Rtowergr + self.Ze)
        self.Istwsum = (self.Istwcond + self.Istwmut)
        self.Uphstw = self.Zphstmut * (self.Istwsum)
    def calculateSpanVoltages(self):
        #positive values are in the Iz direction
        self.Ustwspan = self.Istwsum * self.Zstw - self.Izref * self.Zphstmut
        self.Uphwspan = self.Izref * self.Zphw - self.Istwsum * self.Zphstmut



def exportResultsToLatex(towerList):
    towerid = []
    towerdistance = []
    Istwmut = []
    Istwcond = []
    Istwsum = []
    Ue = []
    Ze = []
    Uphspan = []
    Ustspan = []
    for idx, tower in enumerate(towerList):
        towerid.append(idx)
        towerdistance.append(int(tower.Dtowers*1000))
        Istwmut.append(abs(tower.Istwmut))
        Istwcond.append(abs(tower.Istwcond))
        Istwsum.append(abs(tower.Istwsum))
        Ze.append(abs(tower.Ze))
        Ue.append(abs(tower.Ue))
        Uphspan.append(abs(tower.Uphwspan))
        Ustspan.append(abs(tower.Ustwspan))
    data = {
            "Oszlop": tuple(towerid),
            "Oszlopköz": tuple(towerdistance),
            "Ig indukt.": tuple(Istwmut),
            "Ig kond.": tuple(Istwcond),
            "Ig szumma": tuple(Istwsum),
            "Ue": tuple(Ue),
            "Ze": tuple(Ze),
            "Ustwspan": tuple(Ustspan),
            "Uphwspan": tuple(Uphspan)
        }
    formatteddata = texout.DataFrame(data)
    latextxt = formatteddata.to_latex(
        index = False,
        position = "h",
        float_format="%.2f",
        caption="Távvezeték zárlati áramok nagysága és eloszlása",
        label="tab:tavvezetek_eredmenyek"
        )
    print(latextxt)
    with open("tavvezetekpy.tex", "w") as f:
        f.write(latextxt)

def powerLineCalculation(towerCount):
    powerLine = cPowerLine(towerCount)
    print("Most írom ki a szövegfájlba az adatokat")
    #powerLine.printIntermediateResults()
    #powerLine.printIntermediateResultsOfLine()
    #exportResultsToLatex(powerLine.towerList)
    
def main():
    towerCount = 40  #Number of poles in the line included the faulty one. 
    for towerCount in range(40):
        powerLineCalculation(towerCount+1)

if __name__=="__main__":
    main()
