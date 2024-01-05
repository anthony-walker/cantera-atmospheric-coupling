import math


def KUNKNOWN0(T, M):
    return 3.42857142857143e-33*M + 1.44e-13

def KUNKNOWN4(T, M):
    return (5.77777777777778e-31*M*math.exp(3534/T) + 6.5e-34*M*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*M*math.exp(3534/T) + 1)

def KUNKNOWN7(T, M, H2O):
    return 2.66e-54*H2O*M*math.exp(2200/T) + math.exp(980/T)

def KUNKNOWN8(T, H2O):
    return 3.08e-34*H2O*math.exp(2200/T) + math.exp(600/T)

def KUNKNOWN10(T, M):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*M/T**3.1))**2 + 1))*M/(6.68657815739333e-24*M + 4.0e-12*T**3.1)

def KUNKNOWN15(T, M):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*M*math.exp(-21820/T)))**2 + 1))*M*math.exp(-10650/T)/(4.1e-5*M*math.exp(520/T) + 6.0e+15)

def KUNKNOWN18(T, O2):
    return 1.3e-12*O2*math.exp(-330/T)

def KUNKNOWN20(T, M):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*M*math.exp(-22080/T)/T**3.4))**2 + 1))*M*T**0.1*math.exp(-11000/T)/(607949.833456676*M*math.exp(80/T) + 548352223468120.0*T**3.6)

def KUNKNOWN21(T, O2):
    return 3.3e-39*O2*math.exp(530/T)

def KUNKNOWN23(T, M):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*M/T**1.9))**2 + 1))*M/(9.19166118840122e-28*M*T**0.3 + 2.76761949281346e-10*T**1.6)

def KUNKNOWN25(T, M):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*M/T**2.7))**2 + 1))*M/(6.5211280933241e-25*M*T**0.3 + 1.82662886525688e-10*T**2.4)

def KUNKNOWN26(T, M):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646249e-9*M/T**3.9))**2 + 1))*M*T**0.2/(5.15821743570342e-20*M + 6.07196626492316e-13*T**4.3)

def KUNKNOWN29(T, M):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*M/T**1.26))**2 + 1))*M*T**0.24/(6.75499814951862e-28*M + 5.85084691179052e-12*T**1.74)

def KUNKNOWN31(T, M):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*M/T**4.5))**2 + 1))*M/(4.48947569321853e-19*M + 3.0e-11*T**4.5)

def KUNKNOWN37(T, M):
    return 4.0e-32*M*math.exp(-1000/T)

def KUNKNOWN38(T, O2, N2):
    return 1.5441990796514e-27*N2*O2/T**2.6

def KUNKNOWN39(T, O2):
    return 1.65449901391222e-27*O2**2/T**2.6

def KUNKNOWN40(T, N2):
    return 2.0e-11*N2*math.exp(130/T)

def KUNKNOWN41(T, O2):
    return 3.2e-11*O2*math.exp(67/T)

def KUNKNOWN42(T, H2O):
    return 2.14e-10*H2O

def KUNKNOWN46(T, M):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*M/T**2.6))**2 + 1.0))*M/(6.89374589130092e-25*M + 2.0e-12*T**2.6)

def KUNKNOWN47(T, H2O):
    return 1.2e-15*H2O

def KMT05(T, M):
    return 3.42857142857143e-33*M + 1.44e-13

def KMT06(T, H2O):
    return 1.4e-21*H2O*math.exp(2200/T) + 1

def FC1(T):
    return 0.850000000000000

def K10(T, M):
    return 9.19166118840122e-28*M/T**1.6

def K1I(T):
    return 2.76761949281346e-10/T**0.3

def KR1(T, M):
    return 1.01756082705424e-16*M/T**1.9

def NC1(T):
    return 0.75 - 1.27*math.log10(0.85)

def F1(T, M):
    return 10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*M/T**1.9))**2 + 1))

def KMT01(T, M):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*M/T**1.9))**2 + 1))*M/(9.19166118840122e-28*M*T**0.3 + 2.76761949281346e-10*T**1.6)

def FC2(T):
    return 0.600000000000000

def K20(T, M):
    return 6.75499814951862e-28*M/T**1.5

def K2I(T):
    return 5.85084691179052e-12*T**0.24

def KR2(T, M):
    return 7.47116447301733e-18*M/T**1.26

def NC2(T):
    return 0.75 - 1.27*math.log10(0.6)

def F2(T, M):
    return 10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*M/T**1.26))**2 + 1))

def KMT02(T, M):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*M/T**1.26))**2 + 1))*M*T**0.24/(6.75499814951862e-28*M + 5.85084691179052e-12*T**1.74)

def FC3(T):
    return 0.350000000000000

def K30(T, M):
    return 5.15821743570342e-20*M/T**4.1

def K3I(T):
    return 6.07196626492316e-13*T**0.2

def KR3(T, M):
    return 8.67604494646249e-9*M/T**3.9

def NC3(T):
    return 0.75 - 1.27*math.log10(0.35)

def F3(T, M):
    return 10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646249e-9*M/T**3.9))**2 + 1))

def KMT03(T, M):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646249e-9*M/T**3.9))**2 + 1))*M*T**0.2/(5.15821743570342e-20*M + 6.07196626492316e-13*T**4.3)

def FC4(T):
    return 0.350000000000000

def K40(T, M):
    return 607949.833456676*M*math.exp(-11000/T)/T**3.5

def K4I(T):
    return 548352223468120.0*T**0.1*math.exp(-11080/T)

def KR4(T, M):
    return 3.54310386792477e-10*M*math.exp(-22080/T)/T**3.4

def NC4(T):
    return 0.75 - 1.27*math.log10(0.35)

def F4(T, M):
    return 10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*M*math.exp(-22080/T)/T**3.4))**2 + 1))

def KMT04(T, M):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*M*math.exp(-22080/T)/T**3.4))**2 + 1))*M*T**0.1*math.exp(-11000/T)/(607949.833456676*M*math.exp(80/T) + 548352223468120.0*T**3.6)

def FC7(T):
    return 0.810000000000000

def K70(T, M):
    return 6.5211280933241e-25*M/T**2.4

def K7I(T):
    return 1.82662886525688e-10/T**0.3

def KR7(T, M):
    return 1.09381825613438e-13*M/T**2.7

def NC7(T):
    return 0.75 - 1.27*math.log10(0.81)

def F7(T, M):
    return 10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*M/T**2.7))**2 + 1))

def KMT07(T, M):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*M/T**2.7))**2 + 1))*M/(6.5211280933241e-25*M*T**0.3 + 1.82662886525688e-10*T**2.4)

def FC8(T):
    return 0.410000000000000

def K80(T, M):
    return 4.48947569321853e-19*M/T**4.5

def K8I(T):
    return 3.00000000000000e-11

def KR8(T, M):
    return 1.49649189773951e-8*M/T**4.5

def NC8(T):
    return 0.75 - 1.27*math.log10(0.41)

def F8(T, M):
    return 10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*M/T**4.5))**2 + 1))

def KMT08(T, M):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*M/T**4.5))**2 + 1))*M/(4.48947569321853e-19*M + 3.0e-11*T**4.5)

def FC9(T):
    return 0.400000000000000

def K90(T, M):
    return 6.68657815739333e-24*M/T**3.1

def K9I(T):
    return 4.00000000000000e-12

def KR9(T, M):
    return 1.67164453934833e-12*M/T**3.1

def NC9(T):
    return 0.75 - 1.27*math.log10(0.4)

def F9(T, M):
    return 10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*M/T**3.1))**2 + 1))

def KMT09(T, M):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*M/T**3.1))**2 + 1))*M/(6.68657815739333e-24*M + 4.0e-12*T**3.1)

def FC10(T):
    return 0.400000000000000

def K100(T, M):
    return 4.1e-5*M*math.exp(-10650/T)

def K10I(T):
    return 6.0e+15*math.exp(-11170/T)

def KR10(T, M):
    return 6.83333333333333e-21*M*math.exp(-21820/T)

def NC10(T):
    return 0.75 - 1.27*math.log10(0.4)

def F10(T, M):
    return 10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*M*math.exp(-21820/T)))**2 + 1))

def KMT10(T, M):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*M*math.exp(-21820/T)))**2 + 1))*M*math.exp(-10650/T)/(4.1e-5*M*math.exp(520/T) + 6.0e+15)

def K3(T):
    return 6.5e-34*math.exp(1335/T)

def K4(T):
    return 2.7e-17*math.exp(2199/T)

def K1(T):
    return 2.4e-14*math.exp(460/T)

def K2(T, M):
    return 6.5e-34*M*math.exp(1335/T)/(2.40740740740741e-17*M*math.exp(3534/T) + 1)

def KMT11(T, M):
    return (5.77777777777778e-31*M*math.exp(3534/T) + 6.5e-34*M*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*M*math.exp(3534/T) + 1)

def FC12(T):
    return 0.530000000000000

def K120(T, M):
    return 6.89374589130092e-25*M/T**2.6

def K12I(T):
    return 2.00000000000000e-12

def KR12(T, M):
    return 3.44687294565046e-13*M/T**2.6

def NC12(T):
    return 0.75 - 1.27*math.log10(0.53)

def F12(T, M):
    return 10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*M/T**2.6))**2 + 1.0))

def KMT12(T, M):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*M/T**2.6))**2 + 1.0))*M/(6.89374589130092e-25*M + 2.0e-12*T**2.6)
