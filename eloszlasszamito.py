#!/usr/bin/env python3 
import math, cmath, numpy

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

def main():
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
