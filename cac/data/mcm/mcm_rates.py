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

def KUNKNOWN30(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN31(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN32(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN37(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN38(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN39(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN40(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN45(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN55(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN60(T, RO2):
    return 6.0e-13*RO2

def KUNKNOWN64(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN68(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN75(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN82(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN84(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN88(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN112(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN122(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN124(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN125(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN128(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN133(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN134(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN135(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN137(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN146(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN148(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN152(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN159(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN167(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN169(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN170(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN181(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN182(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN183(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN199(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN201(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN202(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN207(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN208(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN214(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN230(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN231(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN241(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN242(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN243(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN249(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN250(T, H2O):
    return 2.0e-18*H2O

def KUNKNOWN260(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN261(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN292(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN293(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN294(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN306(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN307(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN308(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN323(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN324(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN325(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN337(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN338(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN344(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN355(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN356(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN362(T, H2O):
    return 2.0e-16*H2O

def KUNKNOWN381(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN382(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN383(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN388(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN394(T, H2O):
    return 2.0e-16*H2O

def KUNKNOWN405(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN406(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN416(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN417(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN424(T, H2O):
    return 2.0e-16*H2O

def KUNKNOWN434(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN435(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN436(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN442(T, H2O):
    return 2.0e-16*H2O

def KUNKNOWN464(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN471(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN472(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN473(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN485(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN486(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN487(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN500(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN501(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN502(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN531(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN532(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN549(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN550(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN551(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN574(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN575(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN576(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN583(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN584(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN585(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN590(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN591(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN592(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN611(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN617(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN623(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN633(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN634(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN635(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN641(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN643(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN644(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN651(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN674(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN675(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN676(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN683(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN684(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN685(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN691(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN692(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN693(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN699(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN700(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN701(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN718(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN719(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN729(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN730(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN731(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN741(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN742(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN761(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN763(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN764(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN768(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN770(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN774(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN775(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN776(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN803(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN805(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN806(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN819(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN820(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN821(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN834(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN835(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN836(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN844(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN849(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN850(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN851(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN863(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN864(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN865(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN875(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN876(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN877(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN886(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN887(T, RO2):
    return 2.12e-12*RO2

def KUNKNOWN888(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN902(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN903(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN904(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN917(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN919(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN920(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN933(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN934(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN935(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN942(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN954(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN955(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN956(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN965(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN966(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN980(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN981(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN991(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN993(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN994(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1004(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1005(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1006(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1011(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1015(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1021(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1022(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1023(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1033(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1035(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1036(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1046(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1047(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1048(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1053(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1057(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1063(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1064(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1065(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1073(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN1081(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1094(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1095(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1096(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1106(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1107(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1108(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1118(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1131(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1132(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1142(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1143(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1157(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1158(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1172(T, RO2):
    return 1.0e-13*RO2

def KUNKNOWN1173(T, RO2):
    return 1.8e-12*RO2

def KUNKNOWN1174(T, RO2):
    return 1.0e-13*RO2

def KUNKNOWN1183(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1185(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1190(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN1198(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1211(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1212(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1213(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1223(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1224(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1225(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1235(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1249(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1251(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1252(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1262(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1263(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1264(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1269(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1277(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1278(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1291(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1292(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1293(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1304(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1306(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1307(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1317(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1318(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1319(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1324(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1332(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1334(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1335(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1345(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1346(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1347(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1352(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1360(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1361(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1373(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1374(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1375(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1383(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1384(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1396(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1397(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1398(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1406(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN1414(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1427(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1428(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1429(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1439(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1440(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1441(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1451(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1469(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1471(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1472(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1482(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1483(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN1484(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1490(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN1496(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1504(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1505(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1518(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1519(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1520(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1530(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1531(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1532(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1544(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1545(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1546(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1557(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1558(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1570(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1571(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1572(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1582(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1583(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1593(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1595(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1596(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1606(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1607(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN1608(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1613(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1618(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1619(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1628(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1634(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1635(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1636(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1648(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1649(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1650(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1660(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1671(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1673(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1674(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1684(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1685(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1686(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1691(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1699(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1701(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1702(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1712(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1713(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1714(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1719(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1727(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1728(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1741(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1743(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1744(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1754(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1755(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN1756(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN1762(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN1768(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1776(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN1777(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN1785(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1791(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1792(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1793(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1801(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN1808(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1809(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1819(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1820(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1821(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1830(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1831(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1832(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1842(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1843(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1844(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1863(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1864(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1874(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1875(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1884(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1886(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1887(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1892(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1902(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1904(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1915(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1917(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1918(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1925(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1934(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1936(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1937(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1944(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1951(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1952(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1953(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1963(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1965(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1966(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1973(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1978(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1979(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1980(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1988(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1989(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1990(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2001(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2002(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2007(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2008(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2009(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2012(T, O2, N2):
    return 9.45699740932608e-39*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN2013(T, O2, N2):
    return 1.65237647042071e-38*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN2016(T, O2, N2):
    return 4.71378659280767e-30*10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*(N2 + O2)/T**3.95))**2 + 1))*(N2 + O2)/(4.10746943954162e-21*T**0.85*(N2 + O2) + 1.14761330843492e-9*T**3.1)

def KUNKNOWN2024(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2025(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2038(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2040(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2041(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2046(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2051(T, O2):
    return 2.4e-14*O2*math.exp(-325/T)

def KUNKNOWN2056(T, RO2):
    return 9.74293590248853e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN2057(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN2058(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN2073(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2075(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2080(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2084(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2086(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2090(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2102(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2104(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2108(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2113(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2115(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2116(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2120(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2125(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2126(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2127(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2135(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2142(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2143(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2156(T, O2, N2):
    return 2.92938288982509e-26*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN2157(T, O2, N2):
    return 4.37723880088807e-27*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN2167(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2168(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2181(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2183(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2184(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2204(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2206(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2208(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2219(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2221(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2223(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2228(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2229(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2243(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN2244(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN2250(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2252(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2253(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2257(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2266(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2268(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2269(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2274(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2279(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2284(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2286(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2297(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2299(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2300(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2310(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2311(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2312(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2326(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2328(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2333(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2338(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2339(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2340(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN2347(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2348(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2349(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2355(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2360(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2361(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2362(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2370(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2372(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2373(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2377(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2384(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2393(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2398(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN2406(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2408(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2409(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2413(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2419(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2421(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2426(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2440(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2441(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2442(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2449(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2454(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2455(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2456(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2464(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2466(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2468(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2477(T, RO2):
    return 8.8e-12*RO2

def KUNKNOWN2492(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN2493(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2494(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2500(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2502(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2507(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2516(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2517(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2518(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2529(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2531(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2536(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2540(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2542(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2547(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2553(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN2573(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2575(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2576(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2588(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2589(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2590(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2595(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2596(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2597(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2602(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN2603(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN2608(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2613(T, RO2):
    return 5.52e-14*RO2

def KUNKNOWN2614(T, RO2):
    return 1.84e-14*RO2

def KUNKNOWN2615(T, RO2):
    return 1.84e-14*RO2

def KUNKNOWN2616(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2621(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2622(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2623(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2637(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN2650(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2652(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2653(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2655(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2660(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2661(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2678(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2679(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2686(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2687(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2701(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2703(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2708(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2709(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2714(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2715(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2716(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2721(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN2722(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN2734(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2736(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2737(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2743(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2744(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2745(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2758(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2764(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2765(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2766(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2777(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2779(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2783(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2790(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2792(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2796(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2798(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2800(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2802(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2804(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2806(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2808(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2810(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2812(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2814(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2816(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2821(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2822(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2830(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2832(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2840(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2845(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2849(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2851(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2856(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2862(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2869(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2871(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2872(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2885(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2886(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2887(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2894(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2905(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2906(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2907(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2919(T, O2):
    return 1.8e-14*O2*math.exp(-260/T)

def KUNKNOWN2924(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN2925(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2926(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2941(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2943(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2951(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2956(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2961(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2973(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2975(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2976(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2981(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2986(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2987(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2988(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2993(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3006(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3008(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3009(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3013(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3021(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3023(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3024(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3032(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3033(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3034(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3039(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3052(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3056(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3057(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN3058(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3067(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN3068(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN3069(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN3084(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3085(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3098(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3099(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3100(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3111(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3113(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3114(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3118(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3129(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN3130(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN3131(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN3142(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3152(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3165(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3176(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3180(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3181(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3182(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3189(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3193(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3194(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN3195(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3204(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3205(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3206(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3225(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3235(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3245(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3253(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3262(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3271(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3281(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3292(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3296(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3297(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN3298(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3306(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3307(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3308(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3328(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3330(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3339(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3340(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3341(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3350(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3360(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3361(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3362(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3372(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3383(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3386(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3387(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3392(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3398(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3399(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3411(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3422(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3423(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3424(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3429(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3435(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3436(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3444(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3455(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3457(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3458(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3470(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3481(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3490(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3492(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3494(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3501(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3503(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3505(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3509(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3511(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3519(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3531(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3536(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3540(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3546(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3548(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3550(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3558(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3560(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3564(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3576(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3577(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3578(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3597(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN3598(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN3599(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN3608(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3610(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3614(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3622(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3624(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3660(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3662(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3664(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3669(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3670(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3681(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3683(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3687(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3695(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3702(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3703(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3704(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3714(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3716(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN3718(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3725(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3729(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3745(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3747(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3749(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3751(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3753(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3755(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3757(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3759(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3761(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3763(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3765(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3767(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3769(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3771(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3773(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3775(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3777(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3779(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3781(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3793(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3794(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3795(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3805(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3806(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3820(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3821(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3822(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3840(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3848(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN3865(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3866(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3867(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3882(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3883(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3884(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3896(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3898(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3899(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3907(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3908(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3909(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3916(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3922(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3923(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3924(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3942(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3944(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3945(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3953(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3954(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3955(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3962(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3969(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3971(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3972(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3980(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3981(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3986(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3992(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3993(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3994(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4008(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4010(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4016(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4025(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4026(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN4027(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4037(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4038(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4039(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4049(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN4065(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4067(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4068(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4073(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4078(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4079(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN4080(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4085(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4092(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN4093(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN4110(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4112(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4113(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4118(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4124(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4125(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN4126(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4131(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4136(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4137(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4152(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4154(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4155(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4173(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4178(T, RO2):
    return 1.5e-12*RO2

def KUNKNOWN4179(T, RO2):
    return 5.0e-13*RO2

def KUNKNOWN4180(T, RO2):
    return 5.0e-13*RO2

def KUNKNOWN4187(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4192(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4193(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4207(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4208(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4225(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4227(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4233(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4240(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4244(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4245(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4246(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4254(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN4265(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4266(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN4267(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4276(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4280(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4281(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4282(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4291(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN4292(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4293(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4309(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4310(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4319(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4320(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4321(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4338(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4340(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4341(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4343(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4347(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4348(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN4349(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4356(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4361(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4362(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4371(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4372(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4373(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4385(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4389(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4390(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4391(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4400(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4401(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN4402(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4415(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4417(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4428(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4429(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN4430(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN4443(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4444(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4454(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN4455(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4456(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN4466(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4478(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4480(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4486(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4495(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN4511(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4512(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4513(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4534(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4536(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4537(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4538(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4542(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN4543(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4544(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4550(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4559(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4563(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4564(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4565(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4577(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN4587(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4589(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4590(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4592(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4601(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4606(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4607(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN4608(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4627(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4628(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4629(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4639(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4641(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4643(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4654(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4656(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4657(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4663(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4668(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN4669(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4670(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4683(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4685(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4686(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4689(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4694(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN4695(T, RO2):
    return 5.9663053894349e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN4696(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN4701(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4706(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4708(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4709(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4717(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN4724(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4725(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4726(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4737(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4738(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4753(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4755(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN4757(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4767(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4768(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN4769(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4782(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN4791(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4792(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4793(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4803(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4804(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4805(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4813(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4815(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4817(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4819(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4821(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4823(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4825(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4827(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4829(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4831(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4833(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4835(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4837(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4839(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4841(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4843(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4845(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4847(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4849(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4851(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4853(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4855(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4857(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN4868(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN4873(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN4881(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN4892(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4893(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4903(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN4904(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN4914(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4915(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4916(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4927(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4928(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4929(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4943(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN4951(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN4952(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4953(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN4964(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN4969(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN4970(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4971(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN4983(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN4985(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN4986(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN4998(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN4999(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5000(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5007(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5014(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5015(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5023(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5029(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN5034(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN5035(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN5036(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN5047(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5049(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5050(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5057(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN5061(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5066(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN5080(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5081(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5082(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5092(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5094(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5102(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5106(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5113(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5114(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5131(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5133(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5134(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5145(T, RO2):
    return 5.2e-13*RO2

def KUNKNOWN5146(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5151(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5160(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5161(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5162(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5174(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5175(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5186(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5188(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5189(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5199(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5200(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5201(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5206(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5213(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5215(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5216(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5226(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5227(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5228(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5233(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5240(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5241(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN5242(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5252(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5253(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN5254(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5262(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5263(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5278(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN5279(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN5280(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN5288(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5297(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5305(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5313(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN5327(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5328(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5329(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5341(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5342(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5343(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5354(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN5367(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5369(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5370(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5377(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5378(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5387(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5398(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5399(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5418(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5420(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5422(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5430(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5431(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5444(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5446(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5448(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5453(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5465(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5467(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5469(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5479(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5481(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5482(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5489(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5490(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5491(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5500(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5512(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5514(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5515(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5522(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5523(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5532(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5558(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5559(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5560(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5575(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5576(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5577(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5586(T, RO2):
    return 5.28e-11*RO2

def KUNKNOWN5587(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN5588(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN5603(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5604(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5623(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5627(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5629(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5633(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5642(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5644(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5646(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5652(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5654(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5656(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5662(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5664(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5666(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5675(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN5676(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN5686(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5688(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5689(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5699(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5700(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5701(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5706(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5711(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5712(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5724(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN5725(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN5735(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN5736(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5737(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN5753(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5755(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5758(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN5763(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN5768(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5774(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5776(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5777(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5789(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5790(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5791(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5798(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5816(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN5824(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5832(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN5841(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN5852(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5854(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5855(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5865(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5866(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5867(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5872(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5882(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5884(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5885(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5895(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5896(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN5897(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN5902(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5909(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5910(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN5911(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5921(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5922(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN5923(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5931(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5932(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN5933(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN5947(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN5948(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN5962(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5964(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN5965(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN5972(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN5973(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN5978(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN5985(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN5987(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN5994(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN5998(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6007(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6008(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6021(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6022(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6023(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6033(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6035(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6036(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6046(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6047(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN6048(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6053(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6058(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6063(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6064(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN6065(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6076(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6085(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6086(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6087(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6096(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6098(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6106(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN6114(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN6127(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6129(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6130(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6138(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6139(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN6140(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6145(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6155(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6157(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6158(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6165(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6166(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6167(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6172(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6181(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6183(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6184(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6196(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6197(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6204(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6214(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6215(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN6216(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6225(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6227(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6229(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6236(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6238(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6240(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6245(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6252(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6254(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6256(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6264(T, RO2):
    return 5.28e-11*RO2

def KUNKNOWN6265(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN6266(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN6284(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6286(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6288(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6294(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6301(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6302(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6303(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6311(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6317(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6318(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6319(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6332(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6333(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6346(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6347(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6357(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6359(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6370(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6371(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6378(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6387(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6388(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6398(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN6406(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6408(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6409(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6417(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN6422(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6427(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN6435(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN6443(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6445(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6446(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6456(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN6457(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN6458(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN6463(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6470(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6471(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6472(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6480(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN6487(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN6488(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN6501(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6502(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6511(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6513(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6514(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6524(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN6525(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN6526(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN6531(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6538(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6539(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6540(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6548(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6559(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6560(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6573(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6574(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6575(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6585(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6586(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6587(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6597(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6606(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6608(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6609(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6618(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN6619(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN6620(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN6631(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6632(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6646(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN6647(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN6662(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6663(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN6664(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN6672(T, RO2):
    return 5.28e-11*RO2

def KUNKNOWN6673(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN6674(T, RO2):
    return 1.76e-11*RO2

def KUNKNOWN6686(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6688(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6694(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6710(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6714(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN6715(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6716(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6727(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6729(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6736(T, RO2):
    return 2.58109279182287e-13*RO2*math.exp(1105/T)**0.5

def KUNKNOWN6737(T, RO2):
    return 1.10618262506695e-13*RO2*math.exp(1105/T)**0.5

def KUNKNOWN6741(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN6745(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6753(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6754(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN6755(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6759(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6760(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6761(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN6769(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6771(T, RO2):
    return 3.5e-12*RO2

def KUNKNOWN6772(T, RO2):
    return 1.5e-12*RO2

def KUNKNOWN6781(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6785(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN6786(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6787(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN6794(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6799(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN6800(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN6801(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN6809(T, RO2):
    return 5.52e-14*RO2

def KUNKNOWN6810(T, RO2):
    return 1.84e-14*RO2

def KUNKNOWN6811(T, RO2):
    return 1.84e-14*RO2

def KUNKNOWN6816(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6820(T, RO2):
    return 7.89268015315457e-14*RO2*math.exp(1665/T)**0.5

def KUNKNOWN6821(T, RO2):
    return 2.63089338438486e-14*RO2*math.exp(1665/T)**0.5

def KUNKNOWN6822(T, RO2):
    return 2.63089338438486e-14*RO2*math.exp(1665/T)**0.5

def KUNKNOWN6833(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN6838(T, RO2):
    return 1.6787137933549e-13*RO2*math.exp(1235/T)**0.5

def KUNKNOWN6839(T, RO2):
    return 5.59571264451634e-14*RO2*math.exp(1235/T)**0.5

def KUNKNOWN6840(T, RO2):
    return 5.59571264451634e-14*RO2*math.exp(1235/T)**0.5

def KUNKNOWN6849(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN6850(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6851(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN6857(T, RO2):
    return 6.16e-13*RO2

def KUNKNOWN6858(T, RO2):
    return 2.64e-13*RO2

def KUNKNOWN6867(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN6868(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN6896(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN6906(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN6907(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN6918(T, RO2):
    return 1.58745078663875e-12*RO2

def KUNKNOWN6919(T, RO2):
    return 5.29150262212918e-13*RO2

def KUNKNOWN6920(T, RO2):
    return 5.29150262212918e-13*RO2

def KUNKNOWN6935(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6937(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6942(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN6943(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN6962(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6964(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN6965(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN6972(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN6978(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN6979(T, RO2):
    return 2.00798406368178e-12*RO2

def KUNKNOWN6980(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN6981(T):
    return (1.7e-14*math.exp(1743/T) + 8.8e-12)*math.exp(-1320/T)

def KUNKNOWN6992(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN6994(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN6999(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7002(T, O2):
    return 7.2e-14*O2*math.exp(-1080/T)

def KUNKNOWN7003(T):
    return 1.8924e-10*math.exp(780/T)/(math.exp(1160/T) + 498)

def KUNKNOWN7004(T):
    return 3.8e-13*math.exp(1940/T)/(math.exp(1160/T) + 498)

def KUNKNOWN7007(T, O2, N2):
    return 1.89399755807657e-27*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*(N2 + O2)/T**5.5))**2 + 1))*(N2 + O2)/(1.05222086559809e-16*N2 + 1.05222086559809e-16*O2 + 1.8e-11*T**5.5)

def KUNKNOWN7009(T, RO2):
    return 1.47908e-12*RO2*math.exp(-520/T)

def KUNKNOWN7010(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN7011(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN7012(T, O2, N2):
    return 990000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*(N2 + O2)*math.exp(-20250/T)))**2 + 1))*(N2 + O2)*math.exp(-9690/T)/(9.0e-5*(N2 + O2)*math.exp(870/T) + 1.1e+16)

def KUNKNOWN7015(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7021(T, RO2):
    return 1.02878569196893e-12*RO2

def KUNKNOWN7022(T, RO2):
    return 3.42928563989645e-13*RO2

def KUNKNOWN7023(T, RO2):
    return 3.42928563989645e-13*RO2

def KUNKNOWN7037(T, O2):
    return 1.2e-16*O2*math.exp(1580/T)

def KUNKNOWN7042(T, RO2):
    return 2.99332590941915e-12*RO2

def KUNKNOWN7043(T, RO2):
    return 3.74165738677394e-13*RO2

def KUNKNOWN7044(T, RO2):
    return 3.74165738677394e-13*RO2

def KUNKNOWN7053(T, O2):
    return 3.12e-16*O2*math.exp(1580/T)

def KUNKNOWN7056(T, O2):
    return 1.03e-16*O2*math.exp(1580/T)

def KUNKNOWN7063(T):
    return 2.03512165410849e-10/T**0.9

def KUNKNOWN7066(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7067(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7082(T):
    return 2.03512165410849e-10/T**0.9

def KUNKNOWN7085(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7086(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7093(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7097(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7098(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7099(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7111(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN7112(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN7119(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7121(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7122(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7134(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7135(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7136(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7142(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7148(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7149(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7150(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7160(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN7161(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN7177(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7182(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN7183(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN7184(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN7195(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7197(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7198(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7202(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7206(T, RO2):
    return 6.0e-13*RO2

def KUNKNOWN7211(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7215(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7217(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7222(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7228(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7234(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7236(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7244(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7249(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN7254(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN7262(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7270(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7276(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7280(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7282(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7287(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7290(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7296(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7297(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7298(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7309(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7311(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7312(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7316(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7320(T, O2):
    return 3.5e-12*O2

def KUNKNOWN7321(T, O2):
    return 3.0e-12*O2

def KUNKNOWN7329(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN7330(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN7332(T, RO2):
    return 1.92e-12*RO2

def KUNKNOWN7333(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN7334(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN7335(T, O2):
    return 2.0e-12*O2

def KUNKNOWN7336(T, O2):
    return 3.5e-12*O2

def KUNKNOWN7344(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN7345(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN7347(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN7348(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN7349(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN7391(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7393(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7394(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7398(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7399(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7403(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7404(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7405(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7414(T, RO2):
    return 6.44e-13*RO2

def KUNKNOWN7415(T, RO2):
    return 2.76e-13*RO2

def KUNKNOWN7422(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7426(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7427(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7428(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7443(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7444(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7445(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7454(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7456(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7463(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7469(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN7470(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN7475(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7484(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN7491(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7492(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7493(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7505(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7507(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7508(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7512(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7518(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7520(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7524(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN7525(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7526(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN7559(T, O2, N2):
    return 3.42857142857143e-33*N2 + 3.42857142857143e-33*O2 + 1.44e-13

def KUNKNOWN7593(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN7594(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN7595(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN7608(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN7609(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN7610(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN7624(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7626(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7627(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7632(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7639(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7655(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7657(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7662(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7679(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7681(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7693(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7700(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN7704(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7719(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7730(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7731(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7732(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7740(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN7763(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7765(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7773(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN7774(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN7775(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN7789(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7790(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7791(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7801(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7802(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7803(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7816(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7818(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7825(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN7826(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN7834(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7846(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7848(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7849(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7866(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7868(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7869(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7877(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7878(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN7879(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN7889(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7891(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7905(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN7917(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7919(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7920(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7930(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7932(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN7933(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN7945(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7947(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7952(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7956(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN7958(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN7960(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7980(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN8005(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8007(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8008(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8018(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8019(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8020(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8029(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8031(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN8044(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8046(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8047(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8055(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8063(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8073(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8075(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN8080(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8097(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8098(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8099(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8117(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8119(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8120(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8126(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8135(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8137(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8138(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8143(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8153(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8155(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN8159(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8165(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8167(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN8170(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8174(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN8179(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8191(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN8192(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN8205(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8210(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8211(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8212(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8223(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN8232(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8240(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8241(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8242(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8256(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8261(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8262(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8263(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8278(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8279(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8280(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8286(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8291(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8292(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8293(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8306(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8307(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8308(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8320(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN8321(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN8334(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8335(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8336(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8346(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8347(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN8352(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8354(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8355(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8374(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8399(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8400(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8401(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8409(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8410(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8411(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8422(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8423(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8424(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8439(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8440(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN8441(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8462(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8463(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8464(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8478(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8479(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN8491(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8498(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8507(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8512(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN8513(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8514(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN8522(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8527(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN8528(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8529(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8544(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8545(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8546(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8553(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8554(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN8559(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8561(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8562(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8581(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8602(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN8603(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN8617(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8618(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8619(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8633(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8634(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN8640(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8645(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN8646(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8647(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8653(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN8658(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8659(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8660(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8670(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8671(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8672(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8684(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8690(T, O2):
    return 9.5e-39*O2*math.exp(5270/T)/(7.5e-29*O2*math.exp(5610/T) + 1)

def KUNKNOWN8697(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8698(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN8699(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN8713(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8720(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8729(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8740(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8749(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8758(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8767(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8774(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8785(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8794(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8803(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8812(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8821(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN8834(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN8846(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN8847(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN8858(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8859(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8860(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8877(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8878(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8879(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8895(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8897(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN8898(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN8912(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8913(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN8914(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN8923(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8937(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN8947(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN8949(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN8954(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN8963(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN8964(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN8978(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8979(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN8991(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN8992(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN9006(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9007(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9008(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9019(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN9020(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN9027(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9029(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9038(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9039(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9040(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9052(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN9056(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9061(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN9070(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN9079(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN9084(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9085(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9086(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9105(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9106(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9107(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9114(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9115(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9116(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9126(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9128(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9129(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9151(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9153(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9154(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9159(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9170(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9172(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9173(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9178(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9185(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9187(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9188(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9193(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9206(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9208(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9209(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9211(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9225(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9227(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9228(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9233(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9246(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9248(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9249(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9251(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9258(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9260(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9261(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9266(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9281(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9283(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9289(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9294(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9295(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9296(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9302(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9310(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9311(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9312(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9323(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9325(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9335(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9336(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9337(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9342(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9352(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9353(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9354(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9360(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9364(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN9365(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN9366(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN9377(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9378(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9379(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9393(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN9403(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9404(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9405(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9412(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9417(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9418(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9419(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9428(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9434(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9435(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9436(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9441(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9445(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN9446(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN9456(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN9457(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN9480(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9482(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9493(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9495(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9505(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9506(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9507(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9522(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9524(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9529(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN9534(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9545(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9546(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9547(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9556(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9562(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9563(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9564(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9574(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9576(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9588(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9590(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9591(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9597(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9602(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9603(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9604(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9619(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9621(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9622(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9628(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9633(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9634(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9635(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9651(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9652(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9653(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9665(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN9666(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN9676(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9677(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9678(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9690(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9692(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9693(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9698(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9699(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9700(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9708(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9720(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9722(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9723(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9728(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9732(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9734(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9743(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9745(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9756(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9758(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9759(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9769(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9770(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9771(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9777(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9782(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9783(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN9784(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN9799(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9801(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9816(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9818(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9819(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9836(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9838(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9839(T):
    return 22000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN9840(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9857(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9859(T):
    return 8140000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN9860(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9861(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9870(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9872(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9876(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9893(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9895(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9896(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9904(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9905(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN9906(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN9911(T, O2):
    return 5.0e-12*O2

def KUNKNOWN9913(T, O2):
    return 1.6e-11*O2 - 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN9914(T, O2):
    return 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN9921(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9923(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN9924(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN9931(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9933(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN9937(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN9948(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9949(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9950(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9968(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9969(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9970(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN9979(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN9984(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9985(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN9986(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN9995(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9996(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN9997(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN10014(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10015(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10016(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10027(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10028(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10029(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10038(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10039(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10040(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10055(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10056(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10057(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10076(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10077(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10078(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10090(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10091(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10092(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10099(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN10105(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10106(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10107(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10115(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN10123(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN10137(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10139(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10140(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10145(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN10150(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10151(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10152(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10164(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10165(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10166(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10178(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10180(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10181(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10193(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10194(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10195(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10208(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10210(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10211(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10216(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN10221(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10222(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10223(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10234(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10236(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10237(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10247(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN10248(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN10258(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10259(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10260(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10273(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10275(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10276(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10290(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN10291(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN10300(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN10301(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN10309(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN10315(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN10323(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10324(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN10325(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10336(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10338(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN10339(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10347(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN10348(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN10349(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN10353(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN10366(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN10367(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN10368(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN10374(T, O2, N2):
    return (5.77777777777778e-31*(N2 + O2)*math.exp(3534/T) + 6.5e-34*(N2 + O2)*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*(N2 + O2)*math.exp(3534/T) + 1)

def KUNKNOWN10388(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN10390(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN10399(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10400(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN10401(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10414(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN10415(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN10416(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN10433(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN10438(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN10439(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN10440(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN10451(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN10452(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10453(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN10464(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10465(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10466(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10476(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10477(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10478(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10491(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10492(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN10493(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN10510(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10511(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN10512(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN10515(T, O2, N2, H2O):
    return (2.66e-54*H2O*(N2 + O2)*math.exp(1220/T) + 1)*math.exp(980/T)

def KUNKNOWN10516(T, H2O):
    return 3.08e-34*H2O*math.exp(2200/T) + math.exp(600/T)

def KUNKNOWN10936(T, O2, N2):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*(N2 + O2)/T**3.1))**2 + 1))*(N2 + O2)/(6.68657815739333e-24*N2 + 6.68657815739333e-24*O2 + 4.0e-12*T**3.1)

def KUNKNOWN11158(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11160(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11169(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11170(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11171(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN11183(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11185(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11186(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11195(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11197(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11198(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11206(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11208(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11209(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11215(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11219(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11220(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11221(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11232(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11234(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11235(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11240(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11245(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11246(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11247(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11256(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11257(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11258(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11269(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11270(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11271(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11282(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11283(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN11284(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11296(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11300(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11301(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11302(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11309(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11313(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11314(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11315(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11326(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11327(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11328(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11333(T, O2, N2):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*(N2 + O2)*math.exp(-21820/T)))**2 + 1))*(N2 + O2)*math.exp(-10650/T)/(4.1e-5*(N2 + O2)*math.exp(520/T) + 6.0e+15)

def KUNKNOWN11347(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11348(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11349(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11358(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11359(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11360(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11369(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11370(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11371(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11378(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11380(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11384(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11390(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11392(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11400(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11401(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN11402(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11412(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11414(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11415(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11424(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11425(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11426(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11436(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11437(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11438(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11448(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11449(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11450(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11459(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11460(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11461(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11470(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11471(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN11472(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN11483(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11485(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11495(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11497(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11501(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11503(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11510(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11511(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN11512(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN11520(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11522(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11528(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN11535(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11537(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11543(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN11550(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11552(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11558(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN11566(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN11572(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11574(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11575(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11583(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11585(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11586(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11594(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11596(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11597(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11604(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN11612(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11616(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN11617(T, RO2):
    return 1.07559100033423e-13*RO2*math.exp(1365/T)**0.5

def KUNKNOWN11618(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN11625(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11627(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11628(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11645(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11646(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11647(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN11667(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11669(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11677(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11679(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11680(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11690(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11692(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11698(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN11703(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11706(T, O2):
    return 1.3e-12*O2*math.exp(-330/T)

def KUNKNOWN11720(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11721(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN11722(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11730(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11732(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11733(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11735(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11737(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11741(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11742(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11743(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11750(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN11751(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN11768(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN11769(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN11774(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11778(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11779(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11780(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11788(T, O2):
    return 1.5e-14*O2*math.exp(-230/T)

def KUNKNOWN11792(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN11793(T, RO2):
    return 4.87146795124426e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN11794(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN11800(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11801(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11802(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN11810(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN11814(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11815(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN11816(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN11829(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11831(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11832(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11838(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11839(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN11840(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN11849(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11862(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11864(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11868(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11878(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11880(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN11884(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11897(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11899(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11900(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11904(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11911(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11913(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11914(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11919(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11927(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11929(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11930(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11937(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11939(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11940(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11945(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11952(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11957(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN11958(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN11971(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11973(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11974(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11978(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11983(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11985(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11986(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN11990(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN11995(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN11997(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN11998(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12002(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12009(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12010(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN12011(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12021(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12022(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12023(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12036(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12038(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12039(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12046(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12051(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12052(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN12053(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12061(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12063(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12067(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12071(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12073(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12078(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12084(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN12085(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN12102(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN12113(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN12114(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN12122(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12123(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12124(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12140(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12141(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12142(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12152(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12154(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12155(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12164(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12165(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12166(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12182(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12199(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN12206(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12211(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12212(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12213(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12223(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12227(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12228(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12229(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12242(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN12243(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN12252(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN12253(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN12264(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12265(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12276(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN12280(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN12288(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12290(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12291(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12296(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12300(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12301(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12302(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12308(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12318(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12319(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12320(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12334(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN12335(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN12341(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12343(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12344(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12352(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12353(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12362(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12364(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12365(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12371(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12372(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12382(T, RO2):
    return 6.0e-13*RO2

def KUNKNOWN12383(T, RO2):
    return 1.4e-12*RO2

def KUNKNOWN12390(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12391(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12392(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12403(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN12404(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN12405(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN12413(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12415(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12419(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12426(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN12427(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN12428(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN12433(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12443(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN12444(T, RO2):
    return 2.12e-12*RO2

def KUNKNOWN12445(T, RO2):
    return 2.65e-13*RO2

def KUNKNOWN12456(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN12457(T, RO2):
    return 1.92e-12*RO2

def KUNKNOWN12458(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN12477(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN12478(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN12491(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN12492(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN12493(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN12509(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12510(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN12511(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN12533(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12534(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12535(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12546(T, RO2):
    return 4.4e-14*RO2

def KUNKNOWN12547(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN12548(T, RO2):
    return 4.4e-14*RO2

def KUNKNOWN12559(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12560(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12570(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12571(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12581(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12582(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12583(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12587(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN12588(T, H2O):
    return 2.0e-18*H2O

def KUNKNOWN12596(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12597(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12623(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12624(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12634(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12635(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12636(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12640(T, H2O):
    return 1.4e-17*H2O

def KUNKNOWN12641(T, H2O):
    return 2.0e-18*H2O

def KUNKNOWN12653(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12654(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN12655(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN12666(T, RO2):
    return 4.4e-14*RO2

def KUNKNOWN12667(T, RO2):
    return 7.92e-13*RO2

def KUNKNOWN12668(T, RO2):
    return 4.4e-14*RO2

def KUNKNOWN12682(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12684(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12693(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12697(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12698(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12699(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12710(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12711(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12712(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12725(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN12726(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN12733(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12737(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12738(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12739(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12751(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12755(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12756(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12757(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12765(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12766(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12781(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN12782(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN12791(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12793(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12794(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12807(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12808(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12809(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12819(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN12820(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN12834(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12835(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12836(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12843(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12848(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12849(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12850(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12859(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12863(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12864(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12865(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12877(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN12878(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN12888(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12890(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN12891(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN12897(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12901(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12902(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12903(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12914(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12916(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12920(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN12926(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN12928(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN12946(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN12957(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12958(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12959(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12969(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN12970(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN12983(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12984(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN12985(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN12993(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN12997(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN12998(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN12999(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN13010(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN13011(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN13017(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13019(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13020(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13038(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13040(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13041(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13046(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13050(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13052(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13053(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13063(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13064(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN13070(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13076(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13077(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13084(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN13085(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN13086(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN13093(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13094(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN13108(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13109(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13110(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13127(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13129(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13130(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13134(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13135(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13136(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13142(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13160(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13161(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13162(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13174(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13175(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13176(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13185(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13186(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN13194(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13196(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13197(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13199(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13206(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13208(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13209(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13211(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13213(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13216(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN13223(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13224(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13239(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13240(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13253(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13254(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13263(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13264(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13278(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN13293(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN13300(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN13301(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13302(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13307(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13311(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13312(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13313(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13319(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13320(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13340(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13341(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13342(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13357(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13358(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13359(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13374(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13375(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13387(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN13391(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN13401(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13411(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13416(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13417(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13418(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13428(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13429(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13430(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13438(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13439(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13440(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13452(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13454(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13455(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13459(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13470(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN13491(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13492(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13493(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13510(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13511(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13512(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13520(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13522(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13523(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13532(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN13533(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN13534(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN13550(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13562(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13563(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13579(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN13580(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN13587(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN13588(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13602(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN13618(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN13630(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN13631(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN13641(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13642(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13643(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13648(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13649(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13650(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13658(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13659(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13660(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13670(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN13679(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13680(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13681(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13701(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13702(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13711(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13715(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13716(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN13717(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13727(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN13728(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN13732(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13734(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13735(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13740(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13744(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13746(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13747(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13752(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13759(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13764(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13765(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13766(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13774(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13778(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13779(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13780(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13789(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN13797(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN13801(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN13802(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN13803(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN13815(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13816(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13817(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13827(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13828(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13829(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13835(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13836(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13837(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13857(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13858(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13859(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13870(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN13871(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN13879(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13881(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13882(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13884(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13887(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13888(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13889(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN13893(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13900(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13901(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN13902(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN13926(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13928(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13929(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13936(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13943(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13947(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13948(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13949(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13955(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13959(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13960(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN13961(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN13970(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN13972(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN13973(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN13978(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN13982(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13983(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN13984(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN13990(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN13997(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14001(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14002(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14003(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14013(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14024(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14028(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14043(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN14044(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN14045(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN14057(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN14058(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN14059(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN14068(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14069(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14070(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14077(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14090(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14101(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14105(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14114(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14129(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14130(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14131(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14147(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14148(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14163(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14172(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14174(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN14175(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN14180(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14184(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14185(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14186(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14197(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN14205(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14207(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN14208(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN14221(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14222(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14223(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14232(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN14240(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14245(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14246(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14253(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14254(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14261(T, O2, N2):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*(N2 + O2)*math.exp(-22080/T)/T**3.4))**2 + 1))*T**0.1*(N2 + O2)*math.exp(-11000/T)/(548352223468120.0*T**3.6 + 607949.833456676*(N2 + O2)*math.exp(80/T))

def KUNKNOWN14265(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14266(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN14274(T, RO2):
    return 6.7e-16*RO2

def KUNKNOWN14275(T, RO2):
    return 6.03e-15*RO2

def KUNKNOWN14278(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14282(T, RO2):
    return 2.5e-14*RO2

def KUNKNOWN14283(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN14284(T, RO2):
    return 2.5e-14*RO2

def KUNKNOWN14291(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14292(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14299(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14300(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14306(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14307(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14313(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14314(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14320(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN14321(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN14322(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN14343(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14348(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14349(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14350(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN14351(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14356(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14357(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14358(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN14363(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14364(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14365(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14369(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14372(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN14373(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14374(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN14380(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14384(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14385(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14386(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14397(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14398(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14399(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14407(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14411(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14412(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14413(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14421(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14428(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14436(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN14442(T, RO2):
    return 6.7e-15*RO2

def KUNKNOWN14451(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14452(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14459(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14460(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN14471(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14473(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN14474(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN14479(T, O2):
    return 2.6e-14*O2*math.exp(-255/T)

def KUNKNOWN14483(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN14484(T, RO2):
    return 3.88844441904472e-13*RO2

def KUNKNOWN14485(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN14492(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14498(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14510(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14512(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN14513(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN14525(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14526(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN14533(T, O2):
    return 8.9e-14*O2*math.exp(-550/T)

def KUNKNOWN14538(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14539(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14540(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14549(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14551(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN14555(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN14559(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN14567(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14577(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN14578(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN14587(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14596(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14598(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN14604(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14605(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14613(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN14617(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14620(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14626(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14627(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14634(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14640(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14641(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14649(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14650(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14661(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN14663(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN14667(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN14674(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14684(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14690(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14691(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14700(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14701(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14710(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14711(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14720(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14726(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14732(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14740(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14746(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14752(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14758(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14764(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14768(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14776(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14780(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14788(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14792(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14798(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14805(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14813(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14819(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14820(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14833(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14834(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14835(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14844(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14852(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14858(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14859(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14868(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN14873(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN14877(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14878(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14879(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14885(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14886(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14893(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14894(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14900(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN14901(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN14907(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14908(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN14909(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN14917(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14918(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN14919(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN14926(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14932(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14939(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14945(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14946(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14961(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14968(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN14974(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN14975(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN14984(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN14990(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN14993(T, O2):
    return 3.3e-39*O2*math.exp(530/T)

def KUNKNOWN15046(T, O2, N2):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*(N2 + O2)/T**1.9))**2 + 1))*(N2 + O2)/(9.19166118840122e-28*T**0.3*(N2 + O2) + 2.76761949281346e-10*T**1.6)

def KUNKNOWN15060(T, O2, N2):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*(N2 + O2)/T**2.7))**2 + 1))*(N2 + O2)/(6.5211280933241e-25*T**0.3*(N2 + O2) + 1.82662886525688e-10*T**2.4)

def KUNKNOWN15229(T, O2, N2):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646248e-9*(N2 + O2)/T**3.9))**2 + 1))*T**0.2*(N2 + O2)/(5.15821743570342e-20*N2 + 5.15821743570342e-20*O2 + 6.07196626492316e-13*T**4.3)

def KUNKNOWN15231(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15243(T, O2, N2):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*(N2 + O2)/T**1.26))**2 + 1))*T**0.24*(N2 + O2)/(6.75499814951862e-28*N2 + 6.75499814951862e-28*O2 + 5.85084691179052e-12*T**1.74)

def KUNKNOWN15246(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15247(T, O2, N2):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(4.48947569321853e-19*N2 + 4.48947569321853e-19*O2 + 3.0e-11*T**4.5)

def KUNKNOWN15251(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15252(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15254(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15257(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15260(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15261(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15262(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15267(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15268(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15271(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15272(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15273(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15275(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15277(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15278(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15280(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15281(T, O2, N2):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN15572(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN15573(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN15577(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN15581(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN15585(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15589(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15593(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15594(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15605(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15606(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN15607(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15614(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN15615(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15616(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN15617(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15625(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15626(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN15627(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15634(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN15635(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN15641(T, RO2):
    return 1.0e-13*RO2

def KUNKNOWN15642(T, RO2):
    return 1.8e-12*RO2

def KUNKNOWN15643(T, RO2):
    return 1.0e-13*RO2

def KUNKNOWN15652(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN15660(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15667(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15671(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15675(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15676(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15683(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15688(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15693(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15697(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15698(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15705(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15709(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15713(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15717(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15718(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15725(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15729(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15730(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15737(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15756(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN15758(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15759(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15760(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN15762(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN15763(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN15764(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN15765(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN15766(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN15767(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN15768(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN15776(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN15780(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15784(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15788(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15789(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15796(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15800(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15803(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15807(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15811(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15818(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15819(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15826(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15827(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15831(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15838(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15839(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN15843(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15847(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN15851(T, O2, N2):
    return 4.0e-32*(N2 + O2)*math.exp(-1000/T)

def KUNKNOWN15852(T, O2, N2):
    return 1.5441990796514e-27*N2*O2/T**2.6

def KUNKNOWN15853(T, O2):
    return 1.65449901391222e-27*O2**2/T**2.6

def KUNKNOWN15854(T, N2):
    return 2.0e-11*N2*math.exp(130/T)

def KUNKNOWN15855(T, O2):
    return 3.2e-11*O2*math.exp(67/T)

def KUNKNOWN15856(T, H2O):
    return 2.14e-10*H2O

def KUNKNOWN15913(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN15914(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN15915(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN15924(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15925(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN15926(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN15930(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN15947(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN15948(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN15955(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN15956(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN15961(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN15962(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN15963(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN15977(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN15984(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN15985(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN15986(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN15996(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN15997(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16298(T, O2, N2):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*(N2 + O2)/T**2.6))**2 + 1.0))*(N2 + O2)/(6.89374589130092e-25*N2 + 6.89374589130092e-25*O2 + 2.0e-12*T**2.6)

def KUNKNOWN16465(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16470(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16472(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16474(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16478(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16480(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16489(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16490(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16495(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16496(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16498(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16502(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16503(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16506(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16507(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16508(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN16509(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16511(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16516(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16517(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16524(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16525(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN16526(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16530(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16533(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN16534(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16537(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16538(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16539(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16542(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16543(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16546(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN16547(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16551(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16552(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16553(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16555(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16560(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16561(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16562(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16569(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16570(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16578(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16579(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16580(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16586(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16587(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16592(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16593(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16594(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16598(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16602(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16603(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16604(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16606(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16607(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16608(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN16610(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16611(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16612(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16614(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16615(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16616(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16621(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16622(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN16623(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16626(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16628(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16629(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN16630(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16633(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16634(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16635(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16639(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16640(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16641(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN16642(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16645(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN16652(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16659(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16660(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16664(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16665(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16670(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16671(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN16672(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16683(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16687(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16688(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16692(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN16693(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16697(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN16698(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN16709(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16711(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16712(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN16713(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16714(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16717(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN16718(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN16724(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN16729(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16730(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16731(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16735(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN16736(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16737(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16741(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN16742(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16747(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN16748(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN16750(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16751(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16756(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN16757(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN16762(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16764(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16765(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16768(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16772(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16774(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16776(T, O2, N2):
    return 141100000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*(N2 + O2)*math.exp(-25220/T)))**2 + 1))*(N2 + O2)*math.exp(-11280/T)/(0.0017*(N2 + O2)*math.exp(2660/T) + 8.3e+16)

def KUNKNOWN16778(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16779(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16780(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16786(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16787(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16788(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16791(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16794(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN16797(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN16799(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16802(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16803(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16805(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16809(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16810(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN16811(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16814(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN16816(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16821(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16822(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16823(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16827(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16829(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16830(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN16831(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN16832(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN16834(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16835(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN16836(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN16837(T, RO2):
    return 2.4e-14*RO2

def KUNKNOWN16839(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16841(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16843(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16844(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16845(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN16851(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16855(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16856(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN16857(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN16861(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16863(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16865(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16869(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16871(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16873(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16879(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16880(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16885(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16886(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16888(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16892(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16893(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16896(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16897(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16898(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN16899(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN16901(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16906(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN16907(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN16915(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16916(T, RO2):
    return 5.28e-12*RO2

def KUNKNOWN16917(T, RO2):
    return 1.76e-12*RO2

def KUNKNOWN16923(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16926(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN16927(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN16930(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16931(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16943(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16944(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16946(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN16947(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16948(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16949(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN16951(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN16952(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN16955(T, O2):
    return 1.5e-14*O2*math.exp(-200/T)

def KUNKNOWN16956(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16957(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN16958(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN16961(T, H2O):
    return 1.2e-15*H2O

def KUNKNOWN16963(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN16965(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16969(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16973(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16974(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN16978(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16982(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN16984(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN16985(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16986(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN16989(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN16991(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN16993(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN16996(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN16997(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN16998(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN17002(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN17003(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN17006(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17007(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17009(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN17010(T, O2):
    return 3.0e-12*O2

def KUNKNOWN17011(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN17012(T, O2):
    return 3.5e-12*O2

def KUNKNOWN17014(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN17022(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17023(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17024(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17031(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17032(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17040(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17041(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17042(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17048(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17052(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17053(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17054(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN17055(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN17062(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17063(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17066(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN17067(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17068(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN17069(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17071(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17080(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17081(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17082(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17087(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17088(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17089(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17094(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN17097(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN17098(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN17105(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17106(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17109(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN17110(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17111(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN17112(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17114(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17123(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17124(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17129(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17130(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17131(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17132(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17133(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17146(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17148(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN17150(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN17153(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN17154(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN17160(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17161(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17167(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17168(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN17169(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN17172(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17173(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17174(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17175(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17188(T, O2, N2):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17191(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN17192(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN17197(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17198(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17201(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN17202(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN17206(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN17207(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17210(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN17211(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17212(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN17213(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN17215(T, O2, N2):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN17219(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN17220(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN17221(T, RO2):
    return 5.0e-14*RO2
