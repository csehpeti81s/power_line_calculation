#!/usr/bin/env python3 
import math, cmath, numpy

DEBUG=True
DEBUG=False

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
    if DEBUG:
        print("Aramosztás arányok, ha Z1={:.3f}, Z2={:.3f}"\
            .format(Z1, Z2))
        print("Ir1= {:.3f}, Ir2= {:.3f}\n|Ir1|= {:.3f}, |Ir2|= {:.3f}"\
            .format(Ir1, Ir2, abs(Ir1), abs(Ir2)).replace(".", ","))
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
    if DEBUG:
        print("Másodfokú egyenlet megoldása.")
        print("a: {:.3f}; b: {:.3f}; c: {:.3f}".format(a,b,c).replace(".", ","))
        print("Diszkrimináns: {:.3f}; gyöke: {:.3f}".format(D, cmath.sqrt(D)).replace(".", ","))
        print("x1: {:.3f}".format(x1).replace(".", ","))
        print("x2: {:.3f}".format(x2).replace(".", ","))
    return (x1, x2)


def main():
    I_zref = complex(10,0)
    oszlopkozokSzama = 40
    l_ok = 0.3                                                 #Oszlopköz km-ben
    A_g = 150                                               #mm2-ben megadva
    D_fg = 10                                               #m-ben megadva
    R_sz = complex(5,0)
    R_szall = complex(0,0)
    Z_g = l_ok * complex(0.2992,selfReactance(1000, GMRfromA(A_g)))
    Z_fg = complex(0, l_ok * mutReactance(1000, D_fg))

    print("Z_fg={:.4f}".format(Z_fg))

    #Végtelen távvezeték számítása
    if(False):
        print("Végtelen távvezeték eseténi")
        (Z_e1, Z_e2) = quadraticFormula(1, -Z_g, -Z_g*R_sz)
        Z_e = genuineRoot(Z_e1, Z_e2)
        #print("Z_e: {:.3f}".format(Ze).replace(".", ","))
        #Z_e = replusZ(Ze1, R_sz)+Z_g
        print("Z_e: {:.3f}".format(Z_e).replace(".", ","))
        Ig_ratio, void = currentRatio(Z_e, R_sz)
        txt = f'A zárlatos oszl. Z_e={Z_e:.2f} árameloszlás {abs(Ig_ratio):.2f}.'
        print(txt)
        Ig_ratio_mut = Z_fg/Z_g
        txt = f'A zárlatos oszl árameloszlása {abs(Ig_ratio_mut):.2f}.'
        print(txt)
    
    #Véges távvezeték számítása
    if(True):
        for nOszlopkoz in range(oszlopkozokSzama):
            nOszlopkoz = nOszlopkoz + 1
            print("Meghatározott {} oszlopköz esetén".format(nOszlopkoz))
            Z_e = R_szall + Z_g
            U_e = U_fg = I_zref * Z_fg
            I_sz = abs(U_e / R_sz)
            for i in range(nOszlopkoz + 1):
                Ig_ratio, void = currentRatio(Z_e, R_sz)
                txt = f'A {(nOszlopkoz - i)}. oszl. Z_e={Z_e:.2f} árameloszlás {abs(Ig_ratio):.2f}.'
                #print(txt)
                txt = f'A {nOszlopkoz-i}/{nOszlopkoz}. oszl. Z_e={Z_e:.2f} árameloszlás {abs(Ig_ratio):.2f}, '
                txt2 = f'feszültség: {abs(U_e):.2f} kV, Rsz áram: {I_sz:.2f} kA.'

                Z_e = replusZ(Z_e, R_sz) + Z_g
                U_e = U_e * R_sz / (R_sz + Z_g) + U_fg 
                I_sz = abs(U_e / R_sz)
            print(txt,txt2)

if __name__=="__main__":
    main()
