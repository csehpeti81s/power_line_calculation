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

class cTower():
    def __init__(self):
        self.Izref = 10         #reference value of earth fault current in kA 
        self.Dtowers = 0.300     #lengt in kilometers between this and previous poles
        self.Rtowergr = 5        #earth resistance of the current pole in Ω
        self.Rsubstgr = 0.05    #earth resistance of the substation in Ω
        self.Astw = 95 + 50     #cross section of the static wire in mm2
        self.Aphw = 250 + 40    #cross section of the phase wire in mm2
        self.Rphw = 0.1154      #resistivity of the phase wire in Ω/km
        self.Rstw = 0.2992      #resistivity of the static wire in Ω/km
        self.Rgrr = 0         #resistivity of earth return in Ω/km
        self.Dphstw = 10        #distance of ground and static wire in meter
        self.Dgrr = 1000 #depth of earth return current in meter
        self.Zstw = None        #self impedance of static wire ground loop Ω/span
        self.Zphw = None        #sefl ipmedance of phase wire ground loop Ω/span
        self.Zphstmut = None    #mutual impedance of phase and static to ground loops Ω/span 
        self.Ustphw = None
        self.Uphstw = None
        self.Ue = None 
        self.Ze = None
        self.Istwmut = None
        self.Istwcond = None
        self.calcInternalValues()
    def calcInternalValues(self):
        self.Zstw =  self.Dtowers * complex(
            self.Rstw, selfReactance(self.Dgrr, GMRfromA(self.Astw)))  
        self.Zphw = self.Dtowers * complex(
            self.Rphw, selfReactance(self.Dgrr, GMRfromA(self.Aphw))) 
        self.Zphstmut = self.Dtowers * complex(
            self.Rgrr, mutReactance(self.Dgrr, self.Dphstw))
        self.Ustphw = self.Zphstmut * self.Izref
    def calculateBackward(self, next):
        if next is None:
            self.Ze = self.Zstw + self.Rsubstgr + self.Rgrr 
            self.Ue = self.Ustphw
            return
        self.Ze = self.Zstw + self.Rsubstgr + replusZ(self.Rtowergr, next.Ze)
        self.Ue = self.Ustphw + next.Ue * self.Rtowergr / (self.Rtowergr + next.Ze)
    def calculateForward(self, previous):
        currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        if previous is None:
            self.Istwcond = self.Izref * currRatio
            self.Istwmut = self.Ue / self.Rtowergr
            self.Uphstw = self.Zphstmut * (self.Istwcond + self.Istwmut)
            return
        self.Istwcond = previous.Istwcond * currRatio
        currRatio, void = currentRatio(self.Ze, self.Rtowergr)
        self.Istwmut = currRatio * previous.Istwmut + self.Ue / (self.Rtowergr + self.Ze)
        self.Uphstw = self.Zphstmut * (self.Istwcond + self.Istwmut)

def makePowerLine(towerCount):
    towerList = []              #Towers and spans are counted from the faulty ones. 
    for i in range(towerCount):
        tower = cTower()
        towerList.append(tower) 
    return towerList 

def calculateUeZe(towerList):
    nextTower = None
    for tower in reversed(towerList):
        tower.calculateBackward(nextTower)
        nextTower = tower

def calculateIgrw(towerList):
    previousTower = None
    for tower in towerList:
        tower.calculateForward(previousTower)
        previousTower = tower

def main():
    towerCount = 40     #Number of poles in the line included the faulty one. 
    towerList = makePowerLine(towerCount)
    towerList[-1].Dtowers = 0.1  #Fist tower is closer to substation than others. 
    towerList[-1].calcInternalValues()  #Fist tower is closer to substation than others. 
    calculateUeZe(towerList)
    calculateIgrw(towerList)

    #DEBUG part, should be deleted in production
    headertxt = (f'{"#":>3}{"Táv [m]":>8}{"Zshw [Ω]":>13}{"Zphw [Ω]":>13}'
        f'{"Zmutual [Ω]":>13}{"Umutual [V]":>12}{"Ze [Ω]":>7}{"Ue [V]":>7}'
        f'{"Istwmut [kA]":>13}{"Istwcon [kA]":>13}{"Ist [kA]":>9}'
        f'{"Uphstw [V]":>11}')
    print(headertxt)
    for idx, tower in enumerate(towerList):
        towertxt = (f'{idx:3d}' 
            f'{tower.Dtowers*1000:8.0f}{tower.Zstw:13.2f}'
            f'{tower.Zphw:13.2f}' 
            f'{tower.Zphstmut:13.2f}{abs(tower.Ustphw)*1000:12.2f}'
            f'{abs(tower.Ze):7.2f}{abs(tower.Ue):7.2f}'
            f'{tower.Istwmut:13.2f}{tower.Istwcond:13.2f}{abs(tower.Istwmut+tower.Istwcond):9.2f}'
            f'{abs(tower.Uphstw)*1000:11.2f}')
        print(towertxt)


def tst():
    I_zref = complex(10,0)
    oszlopkozokSzama = 40
    l_ok = 0.3                                                 #Oszlopköz km-ben
    A_g = 95 + 50                                               #mm2-ben megadva
    A_f = 250 + 40
    D_fg = 10                                               #m-ben megadva
    R_sz = complex(5,0)
    R_szall = complex(0.00,0)
    Z_g = l_ok * complex(0.2992,selfReactance(1000, GMRfromA(A_g)))
    Z_f = l_ok * complex(0.1154,selfReactance(1000, GMRfromA(A_f)))
    Z_fg = complex(0, l_ok * mutReactance(1000, D_fg)) 
    U_fg = I_zref * Z_fg

    if(DEBUG):
        debhead= (f'{"":*^80}\n{"Zg":>16}{"Zf":>16}{"Zfg":>16}{"Ufg":>16}\n')
        debtxt = (f'{Z_g:16.3f}{Z_f:16.3f}{Z_fg:16.3f}{U_fg:16.3f}\n{"":*^80}')
        print(debhead,debtxt)

    if(True):
        for nOszlopkoz in range(oszlopkozokSzama):
            pass
        if(True): 
            nOszlopkoz = 40
            nOszlopkoz = nOszlopkoz 
            headertxt = (f'{"n":>3}{"N":>4}{"Ze":>16}{"Ue":>16}{"Ig/Iz\'":>7}{"Isz":>6}')
            print(headertxt)
            Z_e = R_szall + Z_g
            U_e = U_fg
            I_sz = abs(U_e / R_sz)
            for i in range(nOszlopkoz):
                Ig_ratio, void = currentRatio(Z_e, R_sz)
                datatxt = ( f'{nOszlopkoz-i-1:3d} {nOszlopkoz-1:3d}{Z_e:>16.2f}'
                    f'{U_e:>16.2f}{abs(Ig_ratio):7.2f}{abs(I_sz):6.2f}')
                Z_e = replusZ(Z_e, R_sz) + Z_g
                U_e = U_e * R_sz / (R_sz + Z_g) + U_fg 
                I_sz = abs(U_e / R_sz)
                print(datatxt)


if __name__=="__main__":
    main()
