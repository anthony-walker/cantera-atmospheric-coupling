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

def KUNKNOWN1(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN19(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN23(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN30(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN38(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN40(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN41(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN57(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN59(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN60(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN65(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN66(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN72(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN87(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN88(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN89(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN97(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN101(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN102(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN103(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN122(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN123(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN124(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN138(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN143(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN144(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN145(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN157(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN158(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN159(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN169(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN170(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN171(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN180(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN181(T, RO2):
    return 2.12e-12*RO2

def KUNKNOWN182(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN194(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN195(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN196(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN209(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN211(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN212(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN225(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN226(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN227(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN234(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN246(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN247(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN248(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN257(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN258(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN273(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN275(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN276(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN283(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN285(T, O2, N2):
    return 9.45699740932608e-39*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN286(T, O2, N2):
    return 1.65237647042071e-38*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN289(T, O2, N2):
    return 4.71378659280767e-30*10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*(N2 + O2)/T**3.95))**2 + 1))*(N2 + O2)/(4.10746943954162e-21*T**0.85*(N2 + O2) + 1.14761330843492e-9*T**3.1)

def KUNKNOWN297(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN298(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN307(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN309(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN310(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN313(T, O2):
    return 2.4e-14*O2*math.exp(-325/T)

def KUNKNOWN318(T, RO2):
    return 9.74293590248853e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN319(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN320(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN335(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN337(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN341(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN350(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN351(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN352(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN360(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN367(T, O2, N2):
    return 2.92938288982509e-26*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN368(T, O2, N2):
    return 4.37723880088807e-27*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN378(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN380(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN381(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN391(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN393(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN395(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN400(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN401(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN408(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN415(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN416(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN417(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN428(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN429(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN430(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN436(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN441(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN442(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN443(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN448(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN453(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN461(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN463(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN468(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN472(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN474(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN476(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN485(T, RO2):
    return 8.8e-12*RO2

def KUNKNOWN506(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN507(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN522(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN524(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN526(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN528(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN530(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN532(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN534(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN536(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN541(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN542(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN559(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN560(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN561(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN575(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN576(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN577(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN588(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN598(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN611(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN623(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN633(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN643(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN651(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN660(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN669(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN679(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN694(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN705(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN708(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN709(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN714(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN720(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN721(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN733(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN744(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN745(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN746(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN751(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN757(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN758(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN766(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN777(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN779(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN780(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN789(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN791(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN793(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN800(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN802(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN804(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN806(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN817(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN818(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN819(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN838(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN839(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN840(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN883(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN885(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN887(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN896(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN897(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN898(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN909(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN911(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN913(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN920(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN925(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN926(T, RO2):
    return 5.9663053894349e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN927(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN935(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN937(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN938(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN946(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN956(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN965(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN970(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN971(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN972(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN979(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN980(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1010(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1012(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1017(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1018(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1028(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1030(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1031(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1039(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN1040(T, RO2):
    return 2.00798406368178e-12*RO2

def KUNKNOWN1041(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN1042(T):
    return (1.7e-14*math.exp(1743/T) + 8.8e-12)*math.exp(-1320/T)

def KUNKNOWN1049(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1051(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1056(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1059(T, O2):
    return 7.2e-14*O2*math.exp(-1080/T)

def KUNKNOWN1060(T):
    return 1.8924e-10*math.exp(780/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1061(T):
    return 3.8e-13*math.exp(1940/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1064(T, O2, N2):
    return 1.89399755807657e-27*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*(N2 + O2)/T**5.5))**2 + 1))*(N2 + O2)/(1.05222086559809e-16*N2 + 1.05222086559809e-16*O2 + 1.8e-11*T**5.5)

def KUNKNOWN1066(T, RO2):
    return 1.47908e-12*RO2*math.exp(-520/T)

def KUNKNOWN1067(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1068(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1069(T, O2, N2):
    return 990000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*(N2 + O2)*math.exp(-20250/T)))**2 + 1))*(N2 + O2)*math.exp(-9690/T)/(9.0e-5*(N2 + O2)*math.exp(870/T) + 1.1e+16)

def KUNKNOWN1076(T, O2):
    return 1.2e-16*O2*math.exp(1580/T)

def KUNKNOWN1081(T, RO2):
    return 2.99332590941915e-12*RO2

def KUNKNOWN1082(T, RO2):
    return 3.74165738677394e-13*RO2

def KUNKNOWN1083(T, RO2):
    return 3.74165738677394e-13*RO2

def KUNKNOWN1092(T, O2):
    return 3.12e-16*O2*math.exp(1580/T)

def KUNKNOWN1095(T, O2):
    return 1.03e-16*O2*math.exp(1580/T)

def KUNKNOWN1102(T):
    return 2.03512165410849e-10/T**0.9

def KUNKNOWN1105(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1106(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1121(T):
    return 2.03512165410849e-10/T**0.9

def KUNKNOWN1124(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1125(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1136(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1138(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1143(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1146(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1152(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1153(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1154(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1163(T, O2):
    return 3.5e-12*O2

def KUNKNOWN1164(T, O2):
    return 3.0e-12*O2

def KUNKNOWN1172(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1173(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1175(T, RO2):
    return 1.92e-12*RO2

def KUNKNOWN1176(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1177(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1178(T, O2):
    return 2.0e-12*O2

def KUNKNOWN1179(T, O2):
    return 3.5e-12*O2

def KUNKNOWN1187(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1188(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1190(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN1191(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN1192(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN1212(T, O2, N2):
    return 3.42857142857143e-33*N2 + 3.42857142857143e-33*O2 + 1.44e-13

def KUNKNOWN1234(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1236(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1243(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1244(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1252(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1261(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1263(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1272(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1274(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1279(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1286(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1287(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1288(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1300(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1302(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1303(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1309(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1318(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1320(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1321(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1326(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1332(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN1333(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1358(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1363(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1364(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1365(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1375(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1380(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1381(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1382(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1388(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1393(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1394(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1395(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1401(T, O2):
    return 9.5e-39*O2*math.exp(5270/T)/(7.5e-29*O2*math.exp(5610/T) + 1)

def KUNKNOWN1408(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1409(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1410(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1424(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN1433(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN1448(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1450(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1451(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1456(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1458(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1462(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN1463(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1464(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1474(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1480(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1481(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1491(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1492(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1516(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1518(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1538(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1540(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1555(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1557(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1558(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1575(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1577(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1578(T):
    return 22000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1579(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1596(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1598(T):
    return 8140000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1599(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1600(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1609(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1611(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1615(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1626(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1627(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1628(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1633(T, O2):
    return 5.0e-12*O2

def KUNKNOWN1635(T, O2):
    return 1.6e-11*O2 - 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN1636(T, O2):
    return 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN1643(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1645(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1646(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1653(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1655(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1659(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1665(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN1673(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1686(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1688(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1689(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1703(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1704(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1713(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1714(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1722(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1728(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1736(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1737(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1738(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1746(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1747(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1748(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1763(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1764(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1765(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1771(T, O2, N2):
    return (5.77777777777778e-31*(N2 + O2)*math.exp(3534/T) + 6.5e-34*(N2 + O2)*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*(N2 + O2)*math.exp(3534/T) + 1)

def KUNKNOWN1780(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1785(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN1786(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1787(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1798(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1799(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1800(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1810(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1811(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1812(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1816(T, O2, N2, H2O):
    return (2.66e-54*H2O*(N2 + O2)*math.exp(1220/T) + 1)*math.exp(980/T)

def KUNKNOWN1817(T, H2O):
    return 3.08e-34*H2O*math.exp(2200/T) + math.exp(600/T)

def KUNKNOWN1951(T, O2, N2):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*(N2 + O2)/T**3.1))**2 + 1))*(N2 + O2)/(6.68657815739333e-24*N2 + 6.68657815739333e-24*O2 + 4.0e-12*T**3.1)

def KUNKNOWN1988(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1990(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1991(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1997(T, O2, N2):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*(N2 + O2)*math.exp(-21820/T)))**2 + 1))*(N2 + O2)*math.exp(-10650/T)/(4.1e-5*(N2 + O2)*math.exp(520/T) + 6.0e+15)

def KUNKNOWN2002(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2004(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2008(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2014(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2015(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2016(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2025(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2027(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2028(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2036(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2038(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2039(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2047(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2051(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN2052(T, RO2):
    return 1.07559100033423e-13*RO2*math.exp(1365/T)**0.5

def KUNKNOWN2053(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN2060(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2062(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2063(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2087(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2089(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2095(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN2100(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2103(T, O2):
    return 1.3e-12*O2*math.exp(-330/T)

def KUNKNOWN2117(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2118(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2119(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2127(T, O2):
    return 1.5e-14*O2*math.exp(-230/T)

def KUNKNOWN2131(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2132(T, RO2):
    return 4.87146795124426e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2133(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2141(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2143(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2147(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2157(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2159(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2163(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2176(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2178(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2179(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2183(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2190(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2192(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2193(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2198(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2206(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2208(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2209(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2216(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2218(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2219(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2224(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2231(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2236(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN2237(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN2250(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2252(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2253(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2257(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2262(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2264(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2265(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2269(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2274(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2276(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2277(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2281(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2288(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2289(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN2290(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2300(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2301(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2302(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2315(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2317(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2318(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2325(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2330(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2331(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN2332(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2340(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2342(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2346(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2350(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2352(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2357(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2363(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN2364(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN2371(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2373(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2374(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2383(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2384(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2385(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2390(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2400(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN2401(T, RO2):
    return 2.12e-12*RO2

def KUNKNOWN2402(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN2413(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN2414(T, RO2):
    return 1.92e-12*RO2

def KUNKNOWN2415(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN2434(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN2435(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN2448(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN2449(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN2450(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN2466(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2467(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN2468(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2476(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2478(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2482(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2489(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2494(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2496(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2497(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2515(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2517(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2518(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2523(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2527(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2529(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2530(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2540(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2541(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2547(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2553(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2554(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2561(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN2562(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN2563(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN2570(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2571(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2585(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2586(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2587(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2604(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2606(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2607(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2611(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2612(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2613(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2619(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2624(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2626(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2627(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2629(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2637(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2650(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2662(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2670(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2675(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2676(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2677(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2687(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2688(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2689(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2697(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2698(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2699(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2710(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2711(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2723(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2735(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2736(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2740(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2742(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2743(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2748(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2752(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2754(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2755(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2760(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2767(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2772(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2773(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2774(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2782(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2786(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2787(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2788(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2797(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN2805(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN2825(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2826(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2827(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2839(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2840(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2841(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2850(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2851(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2852(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2859(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2869(T, O2, N2):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*(N2 + O2)*math.exp(-22080/T)/T**3.4))**2 + 1))*T**0.1*(N2 + O2)*math.exp(-11000/T)/(548352223468120.0*T**3.6 + 607949.833456676*(N2 + O2)*math.exp(80/T))

def KUNKNOWN2873(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2874(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2879(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2882(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN2883(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2884(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2890(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2894(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2895(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2896(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2907(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2908(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2909(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2917(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2921(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2922(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2923(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2931(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2938(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2944(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2945(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2956(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2958(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2959(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2964(T, O2):
    return 2.6e-14*O2*math.exp(-255/T)

def KUNKNOWN2968(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN2969(T, RO2):
    return 3.88844441904472e-13*RO2

def KUNKNOWN2970(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN2977(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2983(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2995(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2997(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2998(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3007(T, O2):
    return 8.9e-14*O2*math.exp(-550/T)

def KUNKNOWN3012(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3013(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3014(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3023(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3029(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3039(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN3040(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN3049(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3058(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3068(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3074(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3075(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3084(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3090(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3095(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3099(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3100(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3101(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3113(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3116(T, O2):
    return 3.3e-39*O2*math.exp(530/T)

def KUNKNOWN3125(T, O2, N2):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*(N2 + O2)/T**1.9))**2 + 1))*(N2 + O2)/(9.19166118840122e-28*T**0.3*(N2 + O2) + 2.76761949281346e-10*T**1.6)

def KUNKNOWN3129(T, O2, N2):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*(N2 + O2)/T**2.7))**2 + 1))*(N2 + O2)/(6.5211280933241e-25*T**0.3*(N2 + O2) + 1.82662886525688e-10*T**2.4)

def KUNKNOWN3147(T, O2, N2):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646248e-9*(N2 + O2)/T**3.9))**2 + 1))*T**0.2*(N2 + O2)/(5.15821743570342e-20*N2 + 5.15821743570342e-20*O2 + 6.07196626492316e-13*T**4.3)

def KUNKNOWN3149(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3153(T, O2, N2):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*(N2 + O2)/T**1.26))**2 + 1))*T**0.24*(N2 + O2)/(6.75499814951862e-28*N2 + 6.75499814951862e-28*O2 + 5.85084691179052e-12*T**1.74)

def KUNKNOWN3155(T, O2, N2):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(4.48947569321853e-19*N2 + 4.48947569321853e-19*O2 + 3.0e-11*T**4.5)

def KUNKNOWN3156(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3157(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3198(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3199(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3203(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3207(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3211(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3215(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3216(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3226(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3230(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3234(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3238(T, O2, N2):
    return 4.0e-32*(N2 + O2)*math.exp(-1000/T)

def KUNKNOWN3239(T, O2, N2):
    return 1.5441990796514e-27*N2*O2/T**2.6

def KUNKNOWN3240(T, O2):
    return 1.65449901391222e-27*O2**2/T**2.6

def KUNKNOWN3241(T, N2):
    return 2.0e-11*N2*math.exp(130/T)

def KUNKNOWN3242(T, O2):
    return 3.2e-11*O2*math.exp(67/T)

def KUNKNOWN3243(T, H2O):
    return 2.14e-10*H2O

def KUNKNOWN3254(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3255(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3256(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3301(T, O2, N2):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*(N2 + O2)/T**2.6))**2 + 1.0))*(N2 + O2)/(6.89374589130092e-25*N2 + 6.89374589130092e-25*O2 + 2.0e-12*T**2.6)

def KUNKNOWN3324(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3327(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3328(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3329(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3331(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3332(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3333(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3339(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3342(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3343(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3347(T, O2, N2):
    return 141100000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*(N2 + O2)*math.exp(-25220/T)))**2 + 1))*(N2 + O2)*math.exp(-11280/T)/(0.0017*(N2 + O2)*math.exp(2660/T) + 8.3e+16)

def KUNKNOWN3350(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3351(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3353(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3355(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3356(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN3357(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN3358(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN3360(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3361(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN3362(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN3363(T, RO2):
    return 2.4e-14*RO2

def KUNKNOWN3366(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3367(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3368(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3371(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3372(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3376(T, O2):
    return 1.5e-14*O2*math.exp(-200/T)

def KUNKNOWN3377(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3378(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3379(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN3382(T, H2O):
    return 1.2e-15*H2O

def KUNKNOWN3383(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN3384(T, O2):
    return 3.0e-12*O2

def KUNKNOWN3385(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN3386(T, O2):
    return 3.5e-12*O2

def KUNKNOWN3393(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3394(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3395(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3401(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3402(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3410(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3411(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3412(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3418(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3422(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3423(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3424(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3425(T, H2O):
    return 1.0e-17*H2O
