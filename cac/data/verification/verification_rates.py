import math

def KRO2NO(T, M):
    return 2.7e-12 * math.exp(360 / T)

def KRO2HO2(T, M):
    return 2.91e-13 * math.exp(1300 / T)

def KAPHO2(T, M):
    return 5.2e-13 * math.exp(980 / T)

def KAPNO(T, M):
    return 7.5e-12 * math.exp(290 / T)

def KRO2NO3(T, M):
    return 2.3e-12

def KNO3AL(T, M):
    return 1.4e-12 * math.exp(-1860 / T)

def KDEC(T, M):
    return 1.0e6

def KROPRIM(T, M):
    return 2.5e-14 * math.exp(-300 / T)

def KROSEC(T, M):
    return 2.5e-14 * math.exp(-300 / T)

def KCH3O2(T, M):
    return 1.03e-13 * math.exp(365 / T)

def K298CH3O2(T, M):
    return 3.5e-13

def K14ISOM1(T, M):
    return 3.00e7 * math.exp(-5300 / T)

def KD0(T, M):
    return 1.10e-05 * M * math.exp(-10100 / T)

def KDI(T, M):
    return 1.90e17 * math.exp(-14100 / T)

def KRD(T, M):
    return KD0(T, M) / KDI(T, M)

def FCD(T, M):
    return 0.30

def NCD(T, M):
    return 0.75 - 1.27 * (math.log10(FCD(T, M)))

def FD(T, M):
    return 10 ** (math.log10(FCD(T, M)) / (1 + (math.log10(KRD(T, M)) / NCD(T, M)) ** 2))

def KBPAN(T, M):
    return (KD0(T, M) * KDI(T, M) * FD(T, M)) / (KD0(T, M) + KDI(T, M))

def KC0(T, M):
    return 3.28e-28 * M * (T / 300) ** -6.87

def KCI(T, M):
    return 1.125e-11 * (T / 300) ** -1.105

def KRC(T, M):
    return KC0(T, M) / KCI(T, M)

def FCC(T, M):
    return 0.30

def NC(T, M):
    return 0.75 - 1.27 * (math.log10(FCC(T, M)))

def FC(T, M):
    return 10 ** (math.log10(FCC(T, M)) / (1 + (math.log10(KRC(T, M)) / NC(T, M)) ** 2))

def KFPAN(T, M):
    return (KC0(T, M) * KCI(T, M) * FC(T, M)) / (KC0(T, M) + KCI(T, M))

def K10(T, M):
    return 1.0e-31 * M * (T / 300) ** -1.6

def K1I(T, M):
    return 5.0e-11 * (T / 300) ** -0.3

def KR1(T, M):
    return K10(T, M) / K1I(T, M)

def FC1(T, M):
    return 0.85

def NC1(T, M):
    return 0.75 - 1.27 * (math.log10(FC1(T, M)))

def F1(T, M):
    return 10 ** (math.log10(FC1(T, M)) / (1 + (math.log10(KR1(T, M)) / NC1(T, M)) ** 2))

def KMT01(T, M):
    return (K10(T, M) * K1I(T, M) * F1(T, M)) / (K10(T, M) + K1I(T, M))

def K20(T, M):
    return 1.3e-31 * M * (T / 300) ** -1.5

def K2I(T, M):
    return 2.3e-11 * (T / 300) ** 0.24

def KR2(T, M):
    return K20(T, M) / K2I(T, M)

def FC2(T, M):
    return 0.6

def NC2(T, M):
    return 0.75 - 1.27 * (math.log10(FC2(T, M)))

def F2(T, M):
    return 10 ** (math.log10(FC2(T, M)) / (1 + (math.log10(KR2(T, M)) / NC2(T, M)) ** 2))

def KMT02(T, M):
    return (K20(T, M) * K2I(T, M) * F2(T, M)) / (K20(T, M) + K2I(T, M))

def K30(T, M):
    return 3.6e-30 * M * (T / 300) ** -4.1

def K3I(T, M):
    return 1.9e-12 * (T / 300) ** 0.2

def KR3(T, M):
    return K30(T, M) / K3I(T, M)

def FC3(T, M):
    return 0.35

def NC3(T, M):
    return 0.75 - 1.27 * (math.log10(FC3(T, M)))

def F3(T, M):
    return 10 ** (math.log10(FC3(T, M)) / (1 + (math.log10(KR3(T, M)) / NC3(T, M)) ** 2))

def KMT03(T, M):
    return (K30(T, M) * K3I(T, M) * F3(T, M)) / (K30(T, M) + K3I(T, M))

def K40(T, M):
    return 1.3e-3 * M * (T / 300) ** -3.5 * math.exp(-11000 / T)

def K4I(T, M):
    return 9.7e14 * (T / 300) ** 0.1 * math.exp(-11080 / T)

def KR4(T, M):
    return K40(T, M) / K4I(T, M)

def FC4(T, M):
    return 0.35

def NC4(T, M):
    return 0.75 - 1.27 * (math.log10(FC4(T, M)))

def F4(T, M):
    return 10 ** (math.log10(FC4(T, M)) / (1 + (math.log10(KR4(T, M)) / NC4(T, M)) ** 2))

def KMT04(T, M):
    return (K40(T, M) * K4I(T, M) * F4(T, M)) / (K40(T, M) + K4I(T, M))

def KMT05(M):
    return 1.44e-13 * (1 + (M / 4.2e19))

def KMT06(T, H2O):
    return 1 + (1.40e-21 * math.exp(2200 / T) * H2O)

def K70(T, M):
    return 7.4e-31 * M * (T / 300) ** -2.4

def K7I(T, M):
    return 3.3e-11 * (T / 300) ** -0.3

def KR7(T, M):
    return K70(T, M) / K7I(T, M)

def FC7(T, M):
    return 0.81

def NC7(T, M):
    return 0.75 - 1.27 * (math.log10(FC7(T, M)))

def F7(T, M):
    return 10 ** (math.log10(FC7(T, M)) / (1 + (math.log10(KR7(T, M)) / NC7(T, M)) ** 2))

def KMT07(T, M):
    return (K70(T, M) * K7I(T, M) * F7(T, M)) / (K70(T, M) + K7I(T, M))

def K80(T, M):
    return 3.2e-30 * M * (T / 300) ** -4.5

def K8I(T, M):
    return 3.0e-11

def KR8(T, M):
    return K80(T, M) / K8I(T, M)

def FC8(T, M):
    return 0.41

def NC8(T, M):
    return 0.75 - 1.27 * (math.log10(FC8(T, M)))

def F8(T, M):
    return 10 ** (math.log10(FC8(T, M)) / (1 + (math.log10(KR8(T, M)) / NC8(T, M)) ** 2))

def KMT08(T, M):
    return (K80(T, M) * K8I(T, M) * F8(T, M)) / (K80(T, M) + K8I(T, M))

def K90(T, M):
    return 1.4e-31 * M * (T / 300) ** -3.1

def K9I(T, M):
    return 4.0e-12

def KR9(T, M):
    return K90(T, M) / K9I(T, M)

def FC9(T, M):
    return 0.4

def NC9(T, M):
    return 0.75 - 1.27 * (math.log10(FC9(T, M)))

def F9(T, M):
    return 10 ** (math.log10(FC9(T, M)) / (1 + (math.log10(KR9(T, M)) / NC9(T, M)) ** 2))

def KMT09(T, M):
    return (K90(T, M) * K9I(T, M) * F9(T, M)) / (K90(T, M) + K9I(T, M))

def K100(T, M):
    return 4.10e-05 * M * math.exp(-10650 / T)

def K10I(T, M):
    return 6.0e15 * math.exp(-11170 / T)

def KR10(T, M):
    return K100(T, M) / K10I(T, M)

def FC10(T, M):
    return 0.4

def NC10(T, M):
    return 0.75 - 1.27 * (math.log10(FC10(T, M)))

def F10(T, M):
    return 10 ** (math.log10(FC10(T, M)) / (1 + (math.log10(KR10(T, M)) / NC10(T, M)) ** 2))

def KMT10(T, M):
    return (K100(T, M) * K10I(T, M) * F10(T, M)) / (K100(T, M) + K10I(T, M))

def K1(T, M):
    return 2.40e-14 * math.exp(460 / T)

def K3(T, M):
    return 6.50e-34 * math.exp(1335 / T)

def K4(T, M):
    return 2.70e-17 * math.exp(2199 / T)

def K2(T, M):
    return (K3(T, M) * M) / (1 + (K3(T, M) * M) / K4(T, M))

def KMT11(T, M):
    return K1(T, M) + K2(T, M)

def K120(T, M):
    return 2.5e-31 * M * (T / 300) ** -2.6

def K12I(T, M):
    return 2.0e-12

def KR12(T, M):
    return K120(T, M) / K12I(T, M)

def FC12(T, M):
    return 0.53

def NC12(T, M):
    return 0.75 - 1.27 * (math.log10(FC12(T, M)))

def F12(T, M):
    return 10 ** (math.log10(FC12(T, M)) / (1 + (math.log10(KR12(T, M)) / NC12(T, M)) ** 2))

def KMT12(T, M):
    return (K120(T, M) * K12I(T, M) * F12(T, M)) / (K120(T, M) + K12I(T, M))

def K130(T, M):
    return 2.5e-30 * M * (T / 300) ** -5.5

def K13I(T, M):
    return 1.8e-11

def KR13(T, M):
    return K130(T, M) / K13I(T, M)

def FC13(T, M):
    return 0.36

def NC13(T, M):
    return 0.75 - 1.27 * (math.log10(FC13(T, M)))

def F13(T, M):
    return 10 ** (math.log10(FC13(T, M)) / (1 + (math.log10(KR13(T, M)) / NC13(T, M)) ** 2))

def KMT13(T, M):
    return (K130(T, M) * K13I(T, M) * F13(T, M)) / (K130(T, M) + K13I(T, M))

def K140(T, M):
    return 9.0e-5 * math.exp(-9690 / T) * M

def K14I(T, M):
    return 1.1e16 * math.exp(-10560 / T)

def KR14(T, M):
    return K140(T, M) / K14I(T, M)

def FC14(T, M):
    return 0.36

def NC14(T, M):
    return 0.75 - 1.27 * (math.log10(FC14(T, M)))

def F14(T, M):
    return 10 ** (math.log10(FC14(T, M)) / (1 + (math.log10(KR14(T, M)) / NC14(T, M)) ** 2))

def KMT14(T, M):
    return (K140(T, M) * K14I(T, M) * F14(T, M)) / (K140(T, M) + K14I(T, M))

def K150(T, M):
    return 8.6e-29 * M * (T / 300) ** -3.1

def K15I(T, M):
    return 9.0e-12 * (T / 300) ** -0.85

def KR15(T, M):
    return K150(T, M) / K15I(T, M)

def FC15(T, M):
    return 0.48

def NC15(T, M):
    return 0.75 - 1.27 * (math.log10(FC15(T, M)))

def F15(T, M):
    return 10 ** (math.log10(FC15(T, M)) / (1 + (math.log10(KR15(T, M)) / NC15(T, M)) ** 2))

def KMT15(T, M):
    return (K150(T, M) * K15I(T, M) * F15(T, M)) / (K150(T, M) + K15I(T, M))

def K160(T, M):
    return 8.0e-27 * M * (T / 300) ** -3.5

def K16I(T, M):
    return 3.0e-11 * (T / 300) ** -1

def KR16(T, M):
    return K160(T, M) / K16I(T, M)

def FC16(T, M):
    return 0.5

def NC16(T, M):
    return 0.75 - 1.27 * (math.log10(FC16(T, M)))

def F16(T, M):
    return 10 ** (math.log10(FC16(T, M)) / (1 + (math.log10(KR16(T, M)) / NC16(T, M)) ** 2))

def KMT16(T, M):
    return (K160(T, M) * K16I(T, M) * F16(T, M)) / (K160(T, M) + K16I(T, M))

def K170(T, M):
    return 5.0e-30 * M * (T / 300) ** -1.5

def K17I(T, M):
    return 1.0e-12

def KR17(T, M):
    return K170(T, M) / K17I(T, M)

def FC17(T, M):
    return 0.17 * math.exp(-51 / T) + math.exp(-T / 204)

def NC17(T, M):
    return 0.75 - 1.27 * (math.log10(FC17(T, M)))

def F17(T, M):
    return 10 ** (math.log10(FC17(T, M)) / (1 + (math.log10(KR17(T, M)) / NC17(T, M)) ** 2))

def KMT17(T, M):
    return (K170(T, M) * K17I(T, M) * F17(T, M)) / (K170(T, M) + K17I(T, M))

def KMT18(T, O2):
    return 9.5e-39 * O2 * math.exp(5270 / T) / (1 + 7.5e-29 * O2 * math.exp(5610 / T))

def KPPN0(T, M):
    return 1.7e-03 * math.exp(-11280 / T) * M

def KPPNI(T, M):
    return 8.3e16 * math.exp(-13940 / T)

def KRPPN(T, M):
    return KPPN0(T, M) / KPPNI(T, M)

def FCPPN(T, M):
    return 0.36

def NCPPN(T, M):
    return 0.75 - 1.27 * (math.log10(FCPPN(T, M)))

def FPPN(T, M):
    return 10 ** (math.log10(FCPPN(T, M)) / (1 + (math.log10(KRPPN(T, M)) / NCPPN(T, M)) ** 2))

def KBPPN(T, M):
    return (KPPN0(T, M) * KPPNI(T, M) * FPPN(T, M)) / (KPPN0(T, M) + KPPNI(T, M))

def KUNKNOWN2(T, O2):
    return 7.2e-14*O2*math.exp(-1080/T)

def KUNKNOWN3(T):
    return 1.8924e-10*math.exp(780/T)/(math.exp(1160/T) + 498)

def KUNKNOWN4(T):
    return 3.8e-13*math.exp(1940/T)/(math.exp(1160/T) + 498)

def KUNKNOWN7(T, O2, N2):
    return 1.89399755807657e-27*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*(N2 + O2)/T**5.5))**2 + 1))*(N2 + O2)/(1.05222086559809e-16*N2 + 1.05222086559809e-16*O2 + 1.8e-11*T**5.5)

def KUNKNOWN9(T, RO2):
    return 1.47908e-12*RO2*math.exp(-520/T)

def KUNKNOWN10(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN11(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN12(T, O2, N2):
    return 990000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*(N2 + O2)*math.exp(-20250/T)))**2 + 1))*(N2 + O2)*math.exp(-9690/T)/(9.0e-5*(N2 + O2)*math.exp(870/T) + 1.1e+16)

def KUNKNOWN19(T, O2, N2):
    return 3.42857142857143e-33*N2 + 3.42857142857143e-33*O2 + 1.44e-13

def KUNKNOWN27(T, O2, N2):
    return (5.77777777777778e-31*(N2 + O2)*math.exp(3534/T) + 6.5e-34*(N2 + O2)*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*(N2 + O2)*math.exp(3534/T) + 1)

def KUNKNOWN30(T, O2, N2, H2O):
    return (2.66e-54*H2O*(N2 + O2)*math.exp(1220/T) + 1)*math.exp(980/T)

def KUNKNOWN31(T, H2O):
    return 3.08e-34*H2O*math.exp(2200/T) + math.exp(600/T)

def KUNKNOWN33(T, O2, N2):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*(N2 + O2)/T**3.1))**2 + 1))*(N2 + O2)/(6.68657815739333e-24*N2 + 6.68657815739333e-24*O2 + 4.0e-12*T**3.1)

def KUNKNOWN38(T, O2, N2):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*(N2 + O2)*math.exp(-21820/T)))**2 + 1))*(N2 + O2)*math.exp(-10650/T)/(4.1e-5*(N2 + O2)*math.exp(520/T) + 6.0e+15)

def KUNKNOWN41(T, O2):
    return 1.3e-12*O2*math.exp(-330/T)

def KUNKNOWN43(T, O2, N2):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*(N2 + O2)*math.exp(-22080/T)/T**3.4))**2 + 1))*T**0.1*(N2 + O2)*math.exp(-11000/T)/(548352223468120.0*T**3.6 + 607949.833456676*(N2 + O2)*math.exp(80/T))

def KUNKNOWN44(T, O2):
    return 3.3e-39*O2*math.exp(530/T)

def KUNKNOWN46(T, O2, N2):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*(N2 + O2)/T**1.9))**2 + 1))*(N2 + O2)/(9.19166118840122e-28*T**0.3*(N2 + O2) + 2.76761949281346e-10*T**1.6)

def KUNKNOWN48(T, O2, N2):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*(N2 + O2)/T**2.7))**2 + 1))*(N2 + O2)/(6.5211280933241e-25*T**0.3*(N2 + O2) + 1.82662886525688e-10*T**2.4)

def KUNKNOWN49(T, O2, N2):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646248e-9*(N2 + O2)/T**3.9))**2 + 1))*T**0.2*(N2 + O2)/(5.15821743570342e-20*N2 + 5.15821743570342e-20*O2 + 6.07196626492316e-13*T**4.3)

def KUNKNOWN52(T, O2, N2):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*(N2 + O2)/T**1.26))**2 + 1))*T**0.24*(N2 + O2)/(6.75499814951862e-28*N2 + 6.75499814951862e-28*O2 + 5.85084691179052e-12*T**1.74)

def KUNKNOWN54(T, O2, N2):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(4.48947569321853e-19*N2 + 4.48947569321853e-19*O2 + 3.0e-11*T**4.5)

def KUNKNOWN60(T, O2, N2):
    return 4.0e-32*(N2 + O2)*math.exp(-1000/T)

def KUNKNOWN61(T, O2, N2):
    return 1.5441990796514e-27*N2*O2/T**2.6

def KUNKNOWN62(T, O2):
    return 1.65449901391222e-27*O2**2/T**2.6

def KUNKNOWN63(T, N2):
    return 2.0e-11*N2*math.exp(130/T)

def KUNKNOWN64(T, O2):
    return 3.2e-11*O2*math.exp(67/T)

def KUNKNOWN65(T, H2O):
    return 2.14e-10*H2O

def KUNKNOWN69(T, O2, N2):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*(N2 + O2)/T**2.6))**2 + 1.0))*(N2 + O2)/(6.89374589130092e-25*N2 + 6.89374589130092e-25*O2 + 2.0e-12*T**2.6)

def KUNKNOWN70(T, H2O):
    return 1.2e-15*H2O
