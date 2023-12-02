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
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN10(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN17(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN19(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN23(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN30(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN38(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN40(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN41(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN57(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN59(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN60(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN65(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN66(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN72(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN87(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN88(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN89(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN97(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN101(T, RO2):
    return 90332111.4*RO2

def KUNKNOWN102(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN103(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN122(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN123(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN124(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN138(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN143(T, RO2):
    return 90332111.4*RO2

def KUNKNOWN144(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN145(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN157(T, RO2):
    return 963542521.6*RO2

def KUNKNOWN158(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN159(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN169(T, RO2):
    return 1397136656.32*RO2

def KUNKNOWN170(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN171(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN180(T, RO2):
    return 159586730.14*RO2

def KUNKNOWN181(T, RO2):
    return 1276693841.12*RO2

def KUNKNOWN182(T, RO2):
    return 159586730.14*RO2

def KUNKNOWN194(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN195(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN196(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN209(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN211(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN212(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN225(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN226(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN227(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN234(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN246(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN247(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN248(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN257(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN258(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN273(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN275(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN276(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN283(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN285(T, O2, N2):
    return 5.6951369565917e-18*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN286(T, O2, N2):
    return 9.95084369338549e-18*10**(math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*(N2 + O2)/T**1.5) - 0.9525*math.log10((math.exp(51/T) + 0.17*math.exp(T/204))*math.exp(-T/204 - 51/T)))**2 + 1.0))*(N2 + O2)/(2.59807621135332e-26*N2 + 2.59807621135332e-26*O2 + 1.0e-12*T**1.5)

def KUNKNOWN289(T, O2, N2):
    return 2.83870863744886e-9*10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*(N2 + O2)/T**3.95))**2 + 1))*(N2 + O2)/(4.10746943954162e-21*T**0.85*(N2 + O2) + 1.14761330843492e-9*T**3.1)

def KUNKNOWN297(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN298(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN307(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN309(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN310(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN313(T, O2):
    return 14453137.824*O2*math.exp(-325/T)

def KUNKNOWN318(T, RO2):
    return 58673331.4204436*RO2*math.exp(365/T)**0.5

def KUNKNOWN319(T, RO2):
    return 19557777.1401479*RO2*math.exp(365/T)**0.5

def KUNKNOWN320(T, RO2):
    return 19557777.1401479*RO2*math.exp(365/T)**0.5

def KUNKNOWN335(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN337(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN341(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN350(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN351(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN352(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN360(T, RO2):
    return 1204428152.0*RO2

def KUNKNOWN367(T, O2, N2):
    return 1.76411561024623e-5*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN368(T, O2, N2):
    return 2.63603481990816e-6*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(3.74122974434877e-18*T*(N2 + O2) + 9.0e-9*T**3.5)

def KUNKNOWN378(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN380(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN381(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN391(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN393(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN395(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN400(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN401(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN408(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN415(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN416(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN417(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN428(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN429(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN430(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN436(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN441(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN442(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN443(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN448(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN453(T, RO2):
    return 782878298.8*RO2

def KUNKNOWN461(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN463(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN468(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN472(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN474(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN476(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN485(T, RO2):
    return 5299483868.8*RO2

def KUNKNOWN506(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN507(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN522(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN524(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN526(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN528(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN530(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN532(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN534(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN536(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN541(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN542(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN559(T, RO2):
    return 303515894.304*RO2

def KUNKNOWN560(T, RO2):
    return 101171964.768*RO2

def KUNKNOWN561(T, RO2):
    return 101171964.768*RO2

def KUNKNOWN575(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN576(T, RO2):
    return 1397136656.32*RO2

def KUNKNOWN577(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN588(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN598(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN611(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN623(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN633(T, RO2):
    return 1204428152.0*RO2

def KUNKNOWN643(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN651(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN660(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN669(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN679(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN694(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN705(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN708(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN709(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN714(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN720(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN721(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN733(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN744(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN745(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN746(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN751(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN757(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN758(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN766(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN777(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN779(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN780(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN789(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN791(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN793(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN800(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN802(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN804(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN806(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN817(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN818(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN819(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN838(T, RO2):
    return 168619941.28*RO2

def KUNKNOWN839(T, RO2):
    return 168619941.28*RO2

def KUNKNOWN840(T, RO2):
    return 505859823.84*RO2

def KUNKNOWN883(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN885(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN887(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN896(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN897(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN898(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN909(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN911(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN913(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN920(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN925(T, RO2):
    return 11976643.6241079*RO2*math.exp(1985/T)**0.5

def KUNKNOWN926(T, RO2):
    return 35929930.8723236*RO2*math.exp(1985/T)**0.5

def KUNKNOWN927(T, RO2):
    return 11976643.6241079*RO2*math.exp(1985/T)**0.5

def KUNKNOWN935(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN937(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN938(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN946(T, RO2):
    return 150553519.0*RO2

def KUNKNOWN956(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN965(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN970(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN971(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN972(T, RO2):
    return 469726979.28*RO2

def KUNKNOWN979(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN980(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1010(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1012(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1017(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1018(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1028(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1030(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1031(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1039(T, RO2):
    return 403078755.844283*RO2

def KUNKNOWN1040(T, RO2):
    return 1209236267.53285*RO2

def KUNKNOWN1041(T, RO2):
    return 403078755.844283*RO2

def KUNKNOWN1042(T):
    return (10237639.292*math.exp(1743/T) + 5299483868.8)*math.exp(-1320/T)

def KUNKNOWN1049(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1051(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1056(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1059(T, O2):
    return 43359413.472*O2*math.exp(-1080/T)

def KUNKNOWN1060(T):
    return 113962991742.24*math.exp(780/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1061(T):
    return 228841348.88*math.exp(1940/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1064(T, O2, N2):
    return 1.14059198938334e-6*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*(N2 + O2)/T**5.5))**2 + 1))*(N2 + O2)/(1.05222086559809e-16*N2 + 1.05222086559809e-16*O2 + 1.8e-11*T**5.5)

def KUNKNOWN1066(T, RO2):
    return 890722795.53008*RO2*math.exp(-520/T)

def KUNKNOWN1067(T, RO2):
    return 62028049.828*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1068(T, RO2):
    return 62028049.828*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1069(T, O2, N2):
    return 5.9619193524e+32*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*(N2 + O2)*math.exp(-20250/T)))**2 + 1))*(N2 + O2)*math.exp(-9690/T)/(9.0e-5*(N2 + O2)*math.exp(870/T) + 1.1e+16)

def KUNKNOWN1080(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1082(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1087(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1092(T, O2):
    return 2107749266.0*O2

def KUNKNOWN1093(T, O2):
    return 1806642228.0*O2

def KUNKNOWN1101(T):
    return 2.45101128932e+30*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1102(T):
    return 2.45101128932e+30*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1104(T, RO2):
    return 1156251025.92*RO2

def KUNKNOWN1105(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN1106(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN1107(T, O2):
    return 1204428152.0*O2

def KUNKNOWN1108(T, O2):
    return 2107749266.0*O2

def KUNKNOWN1116(T):
    return 6.624354836e+30*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1117(T):
    return 6.624354836e+30*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1119(T, RO2):
    return 963542521.6*RO2

def KUNKNOWN1120(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN1121(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN1141(T, O2, N2):
    return 2.06473397485714e-12*N2 + 2.06473397485714e-12*O2 + 86718826.944

def KUNKNOWN1163(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1165(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1172(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1173(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1181(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1190(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1192(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1201(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1203(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1208(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1215(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN1216(T, RO2):
    return 90332111.4*RO2

def KUNKNOWN1217(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN1229(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1231(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1232(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1238(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1247(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1249(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1250(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1255(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1261(T, RO2):
    return 337239882.56*RO2

def KUNKNOWN1262(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN1290(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN1299(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN1314(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1316(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1317(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1322(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1324(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN1328(T, RO2):
    return 216797067.36*RO2

def KUNKNOWN1329(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN1330(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN1339(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1340(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1350(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1351(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1375(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1377(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1397(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1399(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1414(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1416(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1417(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1434(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1436(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1437(T):
    return 1.3248709672e+31*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1438(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1455(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1457(T):
    return 4.90202257864e+30*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1458(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1459(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1468(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1470(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1474(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1485(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1486(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN1487(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1492(T, O2):
    return 3011070380.0*O2

def KUNKNOWN1494(T, O2):
    return 9635425216.0*O2 - 9635425216.0*O2*math.exp(-550/T)

def KUNKNOWN1495(T, O2):
    return 9635425216.0*O2*math.exp(-550/T)

def KUNKNOWN1502(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1504(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1505(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1512(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1514(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1518(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1524(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN1532(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN1545(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1547(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1548(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1562(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN1563(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN1572(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1573(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN1581(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN1587(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1595(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1596(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN1597(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1605(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1606(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN1607(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1622(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1623(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN1624(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1630(T, O2, N2):
    return (3.47945910577778e-10*(N2 + O2)*math.exp(3534/T) + 3.914391494e-13*(N2 + O2)*math.exp(875/T) + 14453137.824)*math.exp(460/T)/(2.40740740740741e-17*(N2 + O2)*math.exp(3534/T) + 1)

def KUNKNOWN1639(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN1644(T, RO2):
    return 216797067.36*RO2

def KUNKNOWN1645(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN1646(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN1657(T, RO2):
    return 469726979.28*RO2

def KUNKNOWN1658(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN1659(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN1669(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1670(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN1671(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1675(T, O2, N2, H2O):
    return (1.60188944216e-33*H2O*(N2 + O2)*math.exp(1220/T) + 6.02214076e+20)*math.exp(980/T)

def KUNKNOWN1676(T, H2O):
    return (1.85481935408e-13*H2O*math.exp(1600/T) + 6.02214076e+20)*math.exp(600/T)

def KUNKNOWN1808(T, O2, N2):
    return 1.61070059466256e-14*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*(N2 + O2)/T**3.1))**2 + 1))*(N2 + O2)/(6.68657815739333e-24*N2 + 6.68657815739333e-24*O2 + 4.0e-12*T**3.1)

def KUNKNOWN1845(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1847(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1848(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1854(T, O2, N2):
    return 1.48144662696e+32*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*(N2 + O2)*math.exp(-21820/T)))**2 + 1))*(N2 + O2)*math.exp(-10650/T)/(4.1e-5*(N2 + O2)*math.exp(520/T) + 6.0e+15)

def KUNKNOWN1859(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1861(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1865(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1871(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1872(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN1873(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN1882(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1884(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1885(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1893(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1895(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1896(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1904(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN1908(T, RO2):
    return 21591201.3473399*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1909(T, RO2):
    return 64773604.0420194*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1910(T, RO2):
    return 21591201.3473399*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1917(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1919(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN1920(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN1942(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1944(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN1950(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN1955(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1958(T, O2):
    return 782878298.8*O2*math.exp(-330/T)

def KUNKNOWN1972(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1973(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN1974(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN1982(T, O2):
    return 9033211.14*O2*math.exp(-230/T)

def KUNKNOWN1986(T, RO2):
    return 97788885.7007396*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN1987(T, RO2):
    return 293366657.102217*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN1988(T, RO2):
    return 97788885.7007396*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN1996(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1998(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN2002(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2012(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2014(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN2018(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2031(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2033(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2034(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2038(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2045(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2047(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2048(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2053(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2061(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2063(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2064(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2071(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2073(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2074(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2079(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2086(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2091(T, RO2):
    return 385417008.64*RO2

def KUNKNOWN2092(T, RO2):
    return 96354252.16*RO2

def KUNKNOWN2105(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2107(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2108(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2112(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2117(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2119(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2120(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2124(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2129(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2131(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2132(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2136(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2143(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2144(T, RO2):
    return 1397136656.32*RO2

def KUNKNOWN2145(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2155(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2156(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2157(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2170(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2172(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2173(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2180(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2185(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2186(T, RO2):
    return 1397136656.32*RO2

def KUNKNOWN2187(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2195(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2197(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN2201(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2205(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2207(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN2212(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2218(T, RO2):
    return 385417008.64*RO2

def KUNKNOWN2219(T, RO2):
    return 96354252.16*RO2

def KUNKNOWN2226(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2228(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2229(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2238(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2239(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN2240(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2245(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2255(T, RO2):
    return 159586730.14*RO2

def KUNKNOWN2256(T, RO2):
    return 1276693841.12*RO2

def KUNKNOWN2257(T, RO2):
    return 159586730.14*RO2

def KUNKNOWN2268(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN2269(T, RO2):
    return 1156251025.92*RO2

def KUNKNOWN2270(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN2289(T, RO2):
    return 385417008.64*RO2

def KUNKNOWN2290(T, RO2):
    return 96354252.16*RO2

def KUNKNOWN2303(T, RO2):
    return 963542521.6*RO2

def KUNKNOWN2304(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN2305(T, RO2):
    return 120442815.2*RO2

def KUNKNOWN2321(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2322(T, RO2):
    return 1397136656.32*RO2

def KUNKNOWN2323(T, RO2):
    return 174642082.04*RO2

def KUNKNOWN2331(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2333(T, RO2):
    return 6022140760.0*RO2

def KUNKNOWN2337(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2344(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2349(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2351(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2352(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2370(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2372(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2373(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2378(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2382(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2384(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2385(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2395(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2396(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN2402(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2408(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN2409(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN2416(T, RO2):
    return 168619941.28*RO2

def KUNKNOWN2417(T, RO2):
    return 168619941.28*RO2

def KUNKNOWN2418(T, RO2):
    return 505859823.84*RO2

def KUNKNOWN2425(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2426(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN2440(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2441(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2442(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2459(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2461(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2462(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2466(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2467(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2468(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2474(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2479(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2481(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2482(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2484(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2492(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN2505(T, RO2):
    return 1204428152.0*RO2

def KUNKNOWN2517(T, RO2):
    return 1204428152.0*RO2

def KUNKNOWN2525(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN2530(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2531(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN2532(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2542(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2543(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2544(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2552(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2553(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN2554(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2565(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN2566(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2578(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2590(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN2591(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN2595(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2597(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2598(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2603(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2607(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2609(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2610(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2615(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2624(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2632(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2650(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN2651(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2652(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2664(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2665(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN2666(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN2675(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2676(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2677(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2684(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2694(T, O2, N2):
    return 2.00760493699448e+41*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*(N2 + O2)*math.exp(-22080/T)/T**3.4))**2 + 1))*T**0.1*(N2 + O2)*math.exp(-11000/T)/(548352223468120.0*T**3.6 + 607949.833456676*(N2 + O2)*math.exp(80/T))

def KUNKNOWN2698(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2699(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN2704(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN2707(T, RO2):
    return 90332111.4*RO2

def KUNKNOWN2708(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN2709(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN2715(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN2719(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2720(T, RO2):
    return 469726979.28*RO2

def KUNKNOWN2721(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2732(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2733(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2734(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2742(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN2746(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2747(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN2748(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN2756(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN2763(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN2769(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2770(T, H2O):
    return 6022.14076*H2O

def KUNKNOWN2781(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2783(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2784(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2789(T, O2):
    return 15657565.976*O2*math.exp(-255/T)

def KUNKNOWN2793(T, RO2):
    return 78055865.4297456*RO2

def KUNKNOWN2794(T, RO2):
    return 234167596.289237*RO2

def KUNKNOWN2795(T, RO2):
    return 78055865.4297456*RO2

def KUNKNOWN2802(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2808(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2820(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2822(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN2823(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN2832(T, O2):
    return 53597052.764*O2*math.exp(-550/T)

def KUNKNOWN2837(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2838(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2839(T, RO2):
    return 469726979.28*RO2

def KUNKNOWN2848(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN2854(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN2864(T, RO2):
    return 385417008.64*RO2

def KUNKNOWN2865(T, RO2):
    return 96354252.16*RO2

def KUNKNOWN2874(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN2883(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2893(T, RO2):
    return 150553519.0*RO2

def KUNKNOWN2899(T, RO2):
    return 337239882.56*RO2

def KUNKNOWN2900(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN2909(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2915(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2920(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN2924(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2925(T, RO2):
    return 156575659.76*RO2

def KUNKNOWN2926(T, RO2):
    return 469726979.28*RO2

def KUNKNOWN2938(T, RO2):
    return 481771260.8*RO2

def KUNKNOWN2941(T, O2):
    return 1.9873064508e-18*O2*math.exp(530/T)

def KUNKNOWN2950(T, O2, N2):
    return 1.53197363309567e-16*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*(N2 + O2)/T**1.9))**2 + 1))*(N2 + O2)/(9.19166118840122e-28*T**0.3*(N2 + O2) + 2.76761949281346e-10*T**1.6)

def KUNKNOWN2954(T, O2, N2):
    return 7.17338185218158e-14*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*(N2 + O2)/T**2.7))**2 + 1))*(N2 + O2)/(6.5211280933241e-25*T**0.3*(N2 + O2) + 1.82662886525688e-10*T**2.4)

def KUNKNOWN2972(T, O2, N2):
    return 1.88616593706739e-11*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646248e-9*(N2 + O2)/T**3.9))**2 + 1))*T**0.2*(N2 + O2)/(5.15821743570342e-20*N2 + 5.15821743570342e-20*O2 + 6.07196626492316e-13*T**4.3)

def KUNKNOWN2974(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2978(T, O2, N2):
    return 2.38009817676418e-18*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*(N2 + O2)/T**1.26))**2 + 1))*T**0.24*(N2 + O2)/(6.75499814951862e-28*N2 + 6.75499814951862e-28*O2 + 5.85084691179052e-12*T**1.74)

def KUNKNOWN2980(T, O2, N2):
    return 8.11087636894818e-9*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*(N2 + O2)/T**4.5))**2 + 1))*(N2 + O2)/(4.48947569321853e-19*N2 + 4.48947569321853e-19*O2 + 3.0e-11*T**4.5)

def KUNKNOWN2981(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2982(T, O2, N2):
    return 126.421015609926*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*(N2 + O2)/T**7.975))**2 + 1))*(N2 + O2)/(3.41740825082107e-11*T**1.105*(N2 + O2) + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3023(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN3024(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN3028(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3032(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN3036(T, RO2):
    return 150553519.0*RO2

def KUNKNOWN3040(T, RO2):
    return 337239882.56*RO2

def KUNKNOWN3041(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN3051(T, RO2):
    return 529948386.88*RO2

def KUNKNOWN3055(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN3059(T, RO2):
    return 55403694.992*RO2

def KUNKNOWN3063(T, O2, N2):
    return 2.408856304e-11*(N2 + O2)*math.exp(-1000/T)

def KUNKNOWN3064(T, O2, N2):
    return 9.29938421912318e-7*N2*O2/T**2.6

def KUNKNOWN3065(T, O2):
    return 9.96362594906059e-7*O2**2/T**2.6

def KUNKNOWN3066(T, N2):
    return 12044281520.0*N2*math.exp(130/T)

def KUNKNOWN3067(T, O2):
    return 19270850432.0*O2*math.exp(67/T)

def KUNKNOWN3068(T, H2O):
    return 128873812264.0*H2O

def KUNKNOWN3079(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN3080(T, RO2):
    return 722656891.2*RO2

def KUNKNOWN3081(T, RO2):
    return 240885630.4*RO2

def KUNKNOWN3126(T, O2, N2):
    return 8.30302162421714e-16*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*(N2 + O2)/T**2.6))**2 + 1.0))*(N2 + O2)/(6.89374589130092e-25*N2 + 6.89374589130092e-25*O2 + 2.0e-12*T**2.6)

def KUNKNOWN3149(T, RO2):
    return 150553519.0*RO2

def KUNKNOWN3152(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3153(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3154(T, O2, N2):
    return 8.432803706228e+32*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3156(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3157(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN3158(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3164(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3167(T, RO2):
    return 337239882.56*RO2

def KUNKNOWN3168(T, RO2):
    return 144531378.24*RO2

def KUNKNOWN3172(T, O2, N2):
    return 8.49724061236e+34*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*(N2 + O2)*math.exp(-25220/T)))**2 + 1))*(N2 + O2)*math.exp(-11280/T)/(0.0017*(N2 + O2)*math.exp(2660/T) + 8.3e+16)

def KUNKNOWN3175(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN3176(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN3178(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3180(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN3181(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN3182(T, RO2):
    return 216797067.36*RO2

def KUNKNOWN3183(T, RO2):
    return 72265689.12*RO2

def KUNKNOWN3185(T, O2):
    return 15055351.9*O2*math.exp(-300/T)

def KUNKNOWN3186(T, RO2):
    return 4817712.608*RO2

def KUNKNOWN3187(T, RO2):
    return 4817712.608*RO2

def KUNKNOWN3188(T, RO2):
    return 14453137.824*RO2

def KUNKNOWN3191(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3192(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN3193(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3196(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN3197(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN3201(T, O2):
    return 9033211.14*O2*math.exp(-200/T)

def KUNKNOWN3202(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN3203(T, RO2):
    return 30110703.8*RO2

def KUNKNOWN3204(T, RO2):
    return 90332111.4*RO2

def KUNKNOWN3207(T, H2O):
    return 722656.8912*H2O

def KUNKNOWN3208(T, O2):
    return 1505535190.0*O2*math.exp(-480/T)

def KUNKNOWN3209(T, O2):
    return 1806642228.0*O2

def KUNKNOWN3210(T, O2):
    return 1505535190.0*O2*math.exp(-480/T)

def KUNKNOWN3211(T, O2):
    return 2107749266.0*O2

def KUNKNOWN3218(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3219(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN3220(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3226(T, RO2):
    return 4215498532.0*RO2

def KUNKNOWN3227(T, RO2):
    return 1806642228.0*RO2

def KUNKNOWN3235(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3236(T, RO2):
    return 317969032.128*RO2

def KUNKNOWN3237(T, RO2):
    return 105989677.376*RO2

def KUNKNOWN3243(T, O2, N2):
    return 1.25862741884e+33*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*(N2 + O2)*math.exp(-24200/T)))**2 + 1))*(N2 + O2)*math.exp(-10100/T)/(1.1e-5*(N2 + O2)*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3247(T, RO2):
    return 38782586.4944*RO2

def KUNKNOWN3248(T, RO2):
    return 16621108.4976*RO2

def KUNKNOWN3249(T, H2O):
    return 3613.284456*H2O

def KUNKNOWN3250(T, H2O):
    return 6022.14076*H2O
