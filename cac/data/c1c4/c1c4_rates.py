import math


def KUNKNOWN1(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN7(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN9(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN10(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN17(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN19(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN23(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN30(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN38(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

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

def KUNKNOWN72(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

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

def KUNKNOWN196(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN197(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN198(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN211(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN213(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN214(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN227(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN228(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN229(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN236(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN248(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN249(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN250(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN259(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN260(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN275(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN277(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN278(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN285(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN287(T, M):
    return 9.45699740932608e-39*10**(math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*M/T**1.5) - 0.9525*math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T)))**2 + 1.0))*M/(2.59807621135332e-26*M + 1.0e-12*T**1.5)

def KUNKNOWN288(T, M):
    return 1.65237647042071e-38*10**(math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*M/T**1.5) - 0.9525*math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T)))**2 + 1.0))*M/(2.59807621135332e-26*M + 1.0e-12*T**1.5)

def KUNKNOWN291(T, M):
    return 4.71378659280767e-30*10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*M/T**3.95))**2 + 1))*M/(4.10746943954162e-21*M*T**0.85 + 1.14761330843492e-9*T**3.1)

def KUNKNOWN299(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN300(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN309(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN311(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN312(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN315(T, O2):
    return 2.4e-14*O2*math.exp(-325/T)

def KUNKNOWN320(T, RO2):
    return 9.74293590248853e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN321(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN322(T, RO2):
    return 3.24764530082951e-14*RO2*math.exp(365/T)**0.5

def KUNKNOWN337(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN339(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN343(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN352(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN353(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN354(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN362(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN369(T, M):
    return 2.92938288982509e-26*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*M/T**4.5))**2 + 1))*M/(3.74122974434877e-18*M*T + 9.0e-9*T**3.5)

def KUNKNOWN370(T, M):
    return 4.37723880088807e-27*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*M/T**4.5))**2 + 1))*M/(3.74122974434877e-18*M*T + 9.0e-9*T**3.5)

def KUNKNOWN380(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN382(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN383(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN393(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN395(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN397(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN402(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN403(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN410(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN417(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN418(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN419(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN426(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN431(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN432(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN433(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN440(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN441(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN442(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN448(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN453(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN454(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN455(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN460(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN465(T, RO2):
    return 1.3e-12*RO2

def KUNKNOWN473(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN475(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN480(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN484(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN486(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN488(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN497(T, RO2):
    return 8.8e-12*RO2

def KUNKNOWN518(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN519(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN534(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN536(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN538(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN540(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN542(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN544(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN546(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN548(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN550(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN555(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN556(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN573(T, RO2):
    return 5.04e-13*RO2

def KUNKNOWN574(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN575(T, RO2):
    return 1.68e-13*RO2

def KUNKNOWN589(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN590(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN591(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN602(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN612(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN625(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN637(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN647(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN657(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN665(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN674(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN683(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN693(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN708(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN719(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN722(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN723(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN728(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN734(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN735(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN747(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN758(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN759(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN760(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN765(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN771(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN772(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN780(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN791(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN793(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN794(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN803(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN805(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN807(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN814(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN816(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN818(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN820(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN831(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN832(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN833(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN852(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN853(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN854(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN897(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN899(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN901(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN910(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN911(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN912(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN923(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN925(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN927(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN934(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN939(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN940(T, RO2):
    return 5.9663053894349e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN941(T, RO2):
    return 1.98876846314497e-14*RO2*math.exp(1985/T)**0.5

def KUNKNOWN949(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN951(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN952(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN960(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN970(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN982(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN987(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN988(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN989(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN996(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN997(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1021(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1034(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1036(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1041(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1042(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1056(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1058(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1059(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1067(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN1068(T, RO2):
    return 2.00798406368178e-12*RO2

def KUNKNOWN1069(T, RO2):
    return 6.69328021227261e-13*RO2

def KUNKNOWN1070(T):
    return (1.7e-14*math.exp(1743/T) + 8.8e-12)*math.exp(-1320/T)

def KUNKNOWN1077(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1079(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1084(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1087(T, O2):
    return 7.2e-14*O2*math.exp(-1080/T)

def KUNKNOWN1088(T):
    return 1.8924e-10*math.exp(780/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1089(T):
    return 3.8e-13*math.exp(1940/T)/(math.exp(1160/T) + 498)

def KUNKNOWN1092(T, M):
    return 1.89399755807657e-27*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*M/T**5.5))**2 + 1))*M/(1.05222086559809e-16*M + 1.8e-11*T**5.5)

def KUNKNOWN1094(T, RO2):
    return 1.47908e-12*RO2*math.exp(-520/T)

def KUNKNOWN1095(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1096(T, RO2):
    return 1.03e-13*RO2*(math.exp(885/T) - 7.18)*math.exp(-520/T)

def KUNKNOWN1097(T, M):
    return 990000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*M*math.exp(-20250/T)))**2 + 1))*M*math.exp(-9690/T)/(9.0e-5*M*math.exp(870/T) + 1.1e+16)

def KUNKNOWN1108(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1110(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1115(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1120(T, O2):
    return 3.5e-12*O2

def KUNKNOWN1121(T, O2):
    return 3.0e-12*O2

def KUNKNOWN1129(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1130(T):
    return 4070000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1132(T, RO2):
    return 1.92e-12*RO2

def KUNKNOWN1133(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1134(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1135(T, O2):
    return 2.0e-12*O2

def KUNKNOWN1136(T, O2):
    return 3.5e-12*O2

def KUNKNOWN1144(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1145(T):
    return 11000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1147(T, RO2):
    return 1.6e-12*RO2

def KUNKNOWN1148(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN1149(T, RO2):
    return 2.0e-13*RO2

def KUNKNOWN1171(T, M):
    return 3.42857142857143e-33*M + 1.44e-13

def KUNKNOWN1193(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1195(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1202(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1203(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1211(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1220(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1222(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1231(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1233(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1238(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1245(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1246(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN1247(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN1259(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1261(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1262(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1268(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1277(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1279(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1280(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1285(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1291(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN1292(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN1320(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN1329(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN1344(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1346(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1347(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1352(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1354(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1358(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN1359(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1360(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1369(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1370(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1380(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1381(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1405(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1407(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1427(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1429(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1444(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1446(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1447(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1464(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1466(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1467(T):
    return 22000000000.0*math.exp((100000000.0 - 8174*T**2)/T**3)

def KUNKNOWN1468(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1485(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1487(T):
    return 8140000000.0*math.exp((100000000.0 - 8591*T**2)/T**3)

def KUNKNOWN1488(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1489(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1498(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1500(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1504(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1515(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1516(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1517(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1522(T, O2):
    return 5.0e-12*O2

def KUNKNOWN1524(T, O2):
    return 1.6e-11*O2 - 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN1525(T, O2):
    return 1.6e-11*O2*math.exp(-550/T)

def KUNKNOWN1532(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1534(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1535(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1542(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1544(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1548(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1554(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN1562(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN1575(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1577(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1578(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1592(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN1593(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN1602(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1603(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN1611(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN1617(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1625(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1626(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1627(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1635(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1636(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1637(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1652(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1653(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN1654(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN1660(T, M):
    return (5.77777777777778e-31*M*math.exp(3534/T) + 6.5e-34*M*math.exp(875/T) + 2.4e-14)*math.exp(460/T)/(2.40740740740741e-17*M*math.exp(3534/T) + 1)

def KUNKNOWN1669(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1674(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN1675(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1676(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN1687(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN1688(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1689(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN1699(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1700(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1701(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1705(T, M, H2O):
    return 2.66e-54*H2O*M*math.exp(2200/T) + math.exp(980/T)

def KUNKNOWN1706(T, H2O):
    return 3.08e-34*H2O*math.exp(2200/T) + math.exp(600/T)

def KUNKNOWN1857(T, M):
    return 2.67463126295733e-35*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(1.67164453934833e-12*M/T**3.1))**2 + 1))*M/(6.68657815739333e-24*M + 4.0e-12*T**3.1)

def KUNKNOWN1897(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1899(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1900(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1906(T, M):
    return 246000000000.0*10**(math.log10(0.4)/(1.77777777777778*(0.9525*math.log10(0.4) - math.log10(6.83333333333333e-21*M*math.exp(-21820/T)))**2 + 1))*M*math.exp(-10650/T)/(4.1e-5*M*math.exp(520/T) + 6.0e+15)

def KUNKNOWN1911(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1913(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN1917(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN1923(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1924(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN1925(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN1934(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1936(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1937(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1945(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1947(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1948(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1956(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN1960(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1961(T, RO2):
    return 1.07559100033423e-13*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1962(T, RO2):
    return 3.58530333444745e-14*RO2*math.exp(1365/T)**0.5

def KUNKNOWN1969(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1971(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN1972(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1996(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN1998(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN1999(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2005(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2007(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2013(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN2018(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2021(T, O2):
    return 1.3e-12*O2*math.exp(-330/T)

def KUNKNOWN2035(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2036(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2037(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2045(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2047(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2048(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2050(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2052(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2056(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2057(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2058(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2065(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2066(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2083(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2084(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2089(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2093(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2094(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2095(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2103(T, O2):
    return 1.5e-14*O2*math.exp(-230/T)

def KUNKNOWN2107(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2108(T, RO2):
    return 4.87146795124426e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2109(T, RO2):
    return 1.62382265041476e-13*RO2*math.exp(-1835/T)**0.5

def KUNKNOWN2117(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2121(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2122(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2123(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2131(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2133(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2137(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2147(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2149(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2153(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2166(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2168(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2169(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2173(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2180(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2182(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2183(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2188(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2196(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2198(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2199(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2206(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2208(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2209(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2214(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2221(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2226(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN2227(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN2240(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2242(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2243(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2247(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2252(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2254(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2255(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2259(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2264(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2266(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2267(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2271(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2278(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2279(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN2280(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2290(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2291(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2292(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2305(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2307(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2308(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2315(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2320(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2321(T, RO2):
    return 2.32e-12*RO2

def KUNKNOWN2322(T, RO2):
    return 2.9e-13*RO2

def KUNKNOWN2330(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2332(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2336(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2340(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2342(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2347(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2353(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN2354(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN2365(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2367(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2368(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2371(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

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

def KUNKNOWN2390(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

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

def KUNKNOWN2477(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2479(T, RO2):
    return 1.0e-11*RO2

def KUNKNOWN2483(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2490(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2495(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2497(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2498(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2516(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2518(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2519(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2524(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2528(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2530(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2531(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2541(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2542(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2548(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2554(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2555(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2562(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN2563(T, RO2):
    return 2.8e-13*RO2

def KUNKNOWN2564(T, RO2):
    return 8.4e-13*RO2

def KUNKNOWN2571(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2572(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2586(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2587(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2588(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2605(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2607(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2608(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2612(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2613(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2614(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2620(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2625(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2627(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2628(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2630(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2638(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2651(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2663(T, RO2):
    return 2.0e-12*RO2

def KUNKNOWN2671(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2676(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2677(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2678(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2688(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2689(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2690(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2698(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2699(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2700(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2717(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2718(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2730(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2742(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN2743(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN2747(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2749(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2750(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2755(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2759(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2761(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2762(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2767(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2776(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN2784(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN2793(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN2794(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN2802(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2804(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2805(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2807(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2810(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2811(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2812(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2816(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN2832(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2833(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2834(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2846(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2847(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN2848(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN2857(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2858(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2859(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2866(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2876(T, M):
    return 3.33370642933042e+20*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(3.54310386792477e-10*M*math.exp(-22080/T)/T**3.4))**2 + 1))*M*T**0.1*math.exp(-11000/T)/(607949.833456676*M*math.exp(80/T) + 548352223468120.0*T**3.6)

def KUNKNOWN2880(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2881(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2886(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2889(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN2890(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2891(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN2897(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2901(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2902(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN2903(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN2914(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2915(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2916(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2924(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN2928(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2929(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN2930(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN2938(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2945(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN2951(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2952(T, H2O):
    return 1.0e-17*H2O

def KUNKNOWN2963(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN2965(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN2966(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN2971(T, O2):
    return 2.6e-14*O2*math.exp(-255/T)

def KUNKNOWN2975(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN2976(T, RO2):
    return 3.88844441904472e-13*RO2

def KUNKNOWN2977(T, RO2):
    return 1.29614813968157e-13*RO2

def KUNKNOWN2984(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN2990(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3002(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3004(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3005(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3014(T, O2):
    return 8.9e-14*O2*math.exp(-550/T)

def KUNKNOWN3019(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3020(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3021(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3030(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3036(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3046(T, RO2):
    return 6.4e-13*RO2

def KUNKNOWN3047(T, RO2):
    return 1.6e-13*RO2

def KUNKNOWN3056(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3065(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3075(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3081(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3082(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3091(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3097(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3102(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3106(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3107(T, RO2):
    return 2.6e-13*RO2

def KUNKNOWN3108(T, RO2):
    return 7.8e-13*RO2

def KUNKNOWN3120(T, RO2):
    return 8.0e-13*RO2

def KUNKNOWN3123(T, O2):
    return 3.3e-39*O2*math.exp(530/T)

def KUNKNOWN3132(T, M):
    return 2.54390206763561e-37*10**(math.log10(0.85)/(1.77777777777778*(0.9525*math.log10(0.85) - math.log10(1.01756082705424e-16*M/T**1.9))**2 + 1))*M/(9.19166118840122e-28*M*T**0.3 + 2.76761949281346e-10*T**1.6)

def KUNKNOWN3136(T, M):
    return 1.19116808093034e-34*10**(math.log10(0.81)/(1.77777777777778*(0.9525*math.log10(0.81) - math.log10(1.09381825613438e-13*M/T**2.7))**2 + 1))*M/(6.5211280933241e-25*M*T**0.3 + 1.82662886525688e-10*T**2.4)

def KUNKNOWN3159(T, M):
    return 3.13205222567296e-32*10**(math.log10(0.35)/(1.77777777777778*(0.9525*math.log10(0.35) - math.log10(8.67604494646249e-9*M/T**3.9))**2 + 1))*M*T**0.2/(5.15821743570342e-20*M + 6.07196626492316e-13*T**4.3)

def KUNKNOWN3161(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3165(T, M):
    return 3.95224600622617e-39*10**(math.log10(0.6)/(1.77777777777778*(0.9525*math.log10(0.6) - math.log10(7.47116447301733e-18*M/T**1.26))**2 + 1))*M*T**0.24/(6.75499814951862e-28*M + 5.85084691179052e-12*T**1.74)

def KUNKNOWN3167(T, M):
    return 1.34684270796556e-29*10**(math.log10(0.41)/(1.77777777777778*(0.9525*math.log10(0.41) - math.log10(1.49649189773951e-8*M/T**4.5))**2 + 1))*M/(4.48947569321853e-19*M + 3.0e-11*T**4.5)

def KUNKNOWN3168(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3169(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

def KUNKNOWN3214(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3215(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3219(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3223(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3227(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3231(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3232(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3242(T, RO2):
    return 8.8e-13*RO2

def KUNKNOWN3246(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3250(T, RO2):
    return 9.2e-14*RO2

def KUNKNOWN3254(T, M):
    return 4.0e-32*M*math.exp(-1000/T)

def KUNKNOWN3255(T, O2, N2):
    return 1.5441990796514e-27*N2*O2/T**2.6

def KUNKNOWN3256(T, O2):
    return 1.65449901391222e-27*O2**2/T**2.6

def KUNKNOWN3257(T, N2):
    return 2.0e-11*N2*math.exp(130/T)

def KUNKNOWN3258(T, O2):
    return 3.2e-11*O2*math.exp(67/T)

def KUNKNOWN3259(T, H2O):
    return 2.14e-10*H2O

def KUNKNOWN3271(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3272(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3273(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3323(T, M):
    return 1.37874917826018e-36*10**(math.log10(0.53)/(1.77777777777778*(0.9525*math.log10(0.53) - math.log10(3.44687294565046e-13*M/T**2.6))**2 + 1.0))*M/(6.89374589130092e-25*M + 2.0e-12*T**2.6)

def KUNKNOWN3354(T, RO2):
    return 2.5e-13*RO2

def KUNKNOWN3357(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3358(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3359(T, M):
    return 1400300000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3361(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3362(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3363(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3370(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3373(T, RO2):
    return 5.6e-13*RO2

def KUNKNOWN3374(T, RO2):
    return 2.4e-13*RO2

def KUNKNOWN3376(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3379(T, M):
    return 141100000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*M*math.exp(-25220/T)))**2 + 1))*M*math.exp(-11280/T)/(0.0017*M*math.exp(2660/T) + 8.3e+16)

def KUNKNOWN3382(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3383(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3385(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3387(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3388(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN3389(T, RO2):
    return 3.6e-13*RO2

def KUNKNOWN3390(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN3392(T, O2):
    return 2.5e-14*O2*math.exp(-300/T)

def KUNKNOWN3393(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN3394(T, RO2):
    return 8.0e-15*RO2

def KUNKNOWN3395(T, RO2):
    return 2.4e-14*RO2

def KUNKNOWN3397(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3398(T, RO2):
    return 4.0e-13*RO2

def KUNKNOWN3399(T, RO2):
    return 1.2e-12*RO2

def KUNKNOWN3406(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3407(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3408(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3411(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3412(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3416(T, O2):
    return 1.5e-14*O2*math.exp(-200/T)

def KUNKNOWN3417(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3418(T, RO2):
    return 5.0e-14*RO2

def KUNKNOWN3419(T, RO2):
    return 1.5e-13*RO2

def KUNKNOWN3422(T, H2O):
    return 1.2e-15*H2O

def KUNKNOWN3424(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN3425(T, RO2):
    return 4.0e-14*RO2

def KUNKNOWN3426(T, RO2):
    return 1.2e-13*RO2

def KUNKNOWN3430(T, RO2):
    return 2.01e-15*RO2

def KUNKNOWN3431(T, RO2):
    return 4.69e-15*RO2

def KUNKNOWN3433(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN3434(T, O2):
    return 3.0e-12*O2

def KUNKNOWN3435(T, O2):
    return 2.5e-12*O2*math.exp(-480/T)

def KUNKNOWN3436(T, O2):
    return 3.5e-12*O2

def KUNKNOWN3443(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3444(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3445(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3451(T, RO2):
    return 7.0e-12*RO2

def KUNKNOWN3452(T, RO2):
    return 3.0e-12*RO2

def KUNKNOWN3460(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3461(T, RO2):
    return 5.28e-13*RO2

def KUNKNOWN3462(T, RO2):
    return 1.76e-13*RO2

def KUNKNOWN3468(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def KUNKNOWN3472(T, RO2):
    return 6.44e-14*RO2

def KUNKNOWN3473(T, RO2):
    return 2.76e-14*RO2

def KUNKNOWN3474(T, H2O):
    return 6.0e-18*H2O

def KUNKNOWN3475(T, H2O):
    return 1.0e-17*H2O

def K14ISOM1(T):
    return 30000000.0*math.exp(-5300/T)

def K298CH3O2(T):
    return 3.50000000000000e-13

def KAPHO2(T):
    return 5.2e-13*math.exp(980/T)

def KAPNO(T):
    return 7.5e-12*math.exp(290/T)

def KCH3O2(T):
    return 1.03e-13*math.exp(365/T)

def KDEC(T):
    return 1000000.00000000

def KMT05(T, M):
    return 3.42857142857143e-33*M + 1.44e-13

def KMT06(T, H2O):
    return 1.4e-21*H2O*math.exp(2200/T) + 1

def KNO3AL(T):
    return 1.44e-12*math.exp(-1862/T)

def KRO2HO2(T):
    return 2.91e-13*math.exp(1300/T)

def KRO2NO(T):
    return 2.7e-12*math.exp(360/T)

def KRO2NO3(T):
    return 2.30000000000000e-12

def KROPRIM(T):
    return 2.5e-14*math.exp(-300/T)

def KROSEC(T):
    return 2.5e-14*math.exp(-300/T)

def FCD(T):
    return 0.300000000000000

def KD0(T, M):
    return 1.1e-5*M*math.exp(-10100/T)

def KDI(T):
    return 1.9e+17*math.exp(-14100/T)

def KRD(T, M):
    return 5.78947368421053e-23*M*math.exp(-24200/T)

def NCD(T):
    return 0.75 - 1.27*math.log10(0.3)

def FD(T, M):
    return 10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))

def KBPAN(T, M):
    return 2090000000000.0*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(5.78947368421053e-23*M*math.exp(-24200/T)))**2 + 1))*M*math.exp(-10100/T)/(1.1e-5*M*math.exp(4000/T) + 1.9e+17)

def FCPPN(T):
    return 0.360000000000000

def KPPN0(T, M):
    return 0.0017*M*math.exp(-11280/T)

def KPPNI(T):
    return 8.3e+16*math.exp(-13940/T)

def KRPPN(T, M):
    return 2.04819277108434e-20*M*math.exp(-25220/T)

def NCPPN(T):
    return 0.75 - 1.27*math.log10(0.36)

def FPPN(T, M):
    return 10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*M*math.exp(-25220/T)))**2 + 1))

def KBPPN(T, M):
    return 141100000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(2.04819277108434e-20*M*math.exp(-25220/T)))**2 + 1))*M*math.exp(-11280/T)/(0.0017*M*math.exp(2660/T) + 8.3e+16)

def FCC(T):
    return 0.300000000000000

def KC0(T, M):
    return 3.41740825082107e-11*M/T**6.87

def KCI(T):
    return 6.14287260767393e-9/T**1.105

def KRC(T, M):
    return 1658.68274830282*M/T**7.975

def NC(T):
    return 0.75 - 1.27*math.log10(0.3)

def FC(T, M):
    return 10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))

def KFPAN(T, M):
    return 2.09927035332076e-19*10**(math.log10(0.3)/(1.77777777777778*(0.9525*math.log10(0.3) - math.log10(1658.68274830282*M/T**7.975))**2 + 1))*M/(3.41740825082107e-11*M*T**1.105 + 6.14287260767393e-9*T**6.87)

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

def FC13(T):
    return 0.360000000000000

def K130(T, M):
    return 1.05222086559809e-16*M/T**5.5

def K13I(T):
    return 1.80000000000000e-11

def KR13(T, M):
    return 5.84567147554496e-6*M/T**5.5

def NC13(T):
    return 0.75 - 1.27*math.log10(0.36)

def F13(T, M):
    return 10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*M/T**5.5))**2 + 1))

def KMT13(T, M):
    return 1.89399755807657e-27*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(5.84567147554496e-6*M/T**5.5))**2 + 1))*M/(1.05222086559809e-16*M + 1.8e-11*T**5.5)

def FC14(T):
    return 0.360000000000000

def K140(T, M):
    return 9.0e-5*M*math.exp(-9690/T)

def K14I(T):
    return 1.1e+16*math.exp(-10560/T)

def KR14(T, M):
    return 8.18181818181818e-21*M*math.exp(-20250/T)

def NC14(T):
    return 0.75 - 1.27*math.log10(0.36)

def F14(T, M):
    return 10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*M*math.exp(-20250/T)))**2 + 1))

def KMT14(T, M):
    return 990000000000.0*10**(math.log10(0.36)/(1.77777777777778*(0.9525*math.log10(0.36) - math.log10(8.18181818181818e-21*M*math.exp(-20250/T)))**2 + 1))*M*math.exp(-9690/T)/(9.0e-5*M*math.exp(870/T) + 1.1e+16)

def FC15(T):
    return 0.480000000000000

def K150(T, M):
    return 4.10746943954162e-21*M/T**3.1

def K15I(T):
    return 1.14761330843492e-9/T**0.85

def KR15(T, M):
    return 5.81948962075021e-8*M/T**3.95

def NC15(T):
    return 0.75 - 1.27*math.log10(0.48)

def F15(T, M):
    return 10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*M/T**3.95))**2 + 1))

def KMT15(T, M):
    return 4.71378659280767e-30*10**(math.log10(0.48)/(1.77777777777778*(0.9525*math.log10(0.48) - math.log10(5.81948962075021e-8*M/T**3.95))**2 + 1))*M/(4.10746943954162e-21*M*T**0.85 + 1.14761330843492e-9*T**3.1)

def FC16(T):
    return 0.500000000000000

def K160(T, M):
    return 3.74122974434877e-18*M/T**3.5

def K16I(T):
    return 9.0e-9/T

def KR16(T, M):
    return 3.74122974434877e-5*M/T**4.5

def NC16(T):
    return 0.75 - 1.27*math.log10(0.5)

def F16(T, M):
    return 10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*M/T**4.5))**2 + 1))

def KMT16(T, M):
    return 3.3671067699139e-26*10**(math.log10(0.5)/(1.77777777777778*(0.9525*math.log10(0.5) - math.log10(3.74122974434877e-5*M/T**4.5))**2 + 1))*M/(3.74122974434877e-18*M*T + 9.0e-9*T**3.5)

def FC17(T):
    return math.exp(-T/204) + 0.17*math.exp(-51/T)

def K170(T, M):
    return 2.59807621135332e-26*M/T**1.5

def K17I(T):
    return 1.00000000000000e-12

def KR17(T, M):
    return 2.59807621135332e-14*M/T**1.5

def NC17(T):
    return 0.75 - 1.27*math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T))

def F17(T, M):
    return 10**(math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*M/T**1.5) - 0.9525*math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T)))**2 + 1.0))

def KMT17(T, M):
    return 2.59807621135332e-38*10**(math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T))/(1.77777777777778*(math.log10(2.59807621135332e-14*M/T**1.5) - 0.9525*math.log10(math.exp(-T/204) + 0.17*math.exp(-51/T)))**2 + 1.0))*M/(2.59807621135332e-26*M + 1.0e-12*T**1.5)
