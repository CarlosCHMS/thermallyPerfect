import numpy
import matplotlib.pyplot as plt

class gasBasic():

    def cp(self, T):

        if T < 1000:
            a = self.a0

        else:
            a = self.a1

        cp = 0.0        
        for ii in range(0, 8):
            cp += a[ii]*(T**(ii-2))

        cp *= self.R/self.M

        return cp

    def toArray(self):
        self.a0 = numpy.array(self.a0)
        self.a1 = numpy.array(self.a1)

        return None

class gasN2(gasBasic):

    def __init__(self):

        self.R = 8.31446261815324
        self.frac = 0.75530
        self.M = 0.0280160
        self.a0 = [-1.33984200E+01, 1.34280300E+00, 3.45742000E+00, 5.74727600E-04,-3.21711900E-06, 7.50775400E-09, -5.90150500E-12, 1.50979900E-15]
        self.a1 = [5.87702841E+05, -2.23921563E+03, 6.06686971E+00, -6.13957913E-04, 1.49178026E-07, -1.92307130E-11, 1.06193594E-15, 0.00000000E+00]

class gasO2(gasBasic):

    def __init__(self):

        self.R = 8.31446261815324
        self.frac = 0.23140
        self.M = 0.032
        self.a0 = [3.88517500E+01, -2.70630800E+00, 3.56119600E+00, -3.32782400E-04, -1.18148000E-06, 1.10853500E-08, -1.49299400E-11, 5.99553800E-15]
        self.a1 = [-1.05642070E+06, 2.41123849E+03, 1.73474238E+00, 1.31512292E-03,-2.29995151E-07, 2.13144378E-11, -7.87498771E-16, 0.00000000E+00]

class gasAr(gasBasic):

    def __init__(self):

        self.R = 8.31446261815324
        self.frac = 0.01290
        self.M = 0.039944
        self.a0 = [0, 0, 2.5, 0, 0, 0 ,0 , 0]
        self.a1 = [0, 0, 2.5, 0, 0, 0 ,0 , 0]

class gasCO2(gasBasic):

    def __init__(self):

        self.R = 8.31446261815324
        self.frac = 0.00040
        self.M = 0.039944
        self.a0 = [-5.34531900E+03, 2.28785400E+02, 2.56174500E-01, 1.68213800E-02,-2.09949800E-05, 1.40430800E-08, -3.81732100E-12, 0.00000000E+00]
        self.a1 = [1.15460081E+05, -1.78337370E+03, 8.28644239E+00,  -8.98356945E-05, 4.26107946E-09,  -1.81443266E-12,  6.29130739E-16, 0.00000000E+00]


class gas():

    def __init__(self):

        """
        Based on:

        Witte, D. W. and Tatum K. E., Computer Code for Determination Perfect Gas Properties of Thermally Perfect Gas Properties, 1994


        """
        self.g = []
        self.g.append(gasN2())
        self.g.append(gasO2())
        self.g.append(gasAr())
        self.g.append(gasCO2())
    
        self.Rgas = 287.035


    def cp(self, T):

        cp = 0
        for g in self.g:
            cp += g.frac*g.cp(T)

        return cp

    def cv(self, T):

        return self.cp(T) - self.Rgas
    
    def toArray(self):
        for g in self.g:
            g.toArray()
    
        return None

    def calcCoeffAir(self):

        self.toArray()

        self.a0 = 0
        self.a1 = 0
        for g in self.g:
            self.a0 += g.a0*g.frac*g.R/g.M
            self.a1 += g.a1*g.frac*g.R/g.M

        print(self.a0, self.a1)


    def cp2(self, T):

        if T < 1000:
            a = self.a0

        else:
            a = self.a1

        cp = 0.0        
        for ii in range(0, 8):
            cp += a[ii]*(T**(ii-2))

        return cp

    def enthalpy(self, T):

        if T < 1000:
            a = self.a0

        else:
            a = self.a1


        h = 0.0        
        for ii in range(2, 8):
            h += a[ii]*( (T**(ii-1)) - (300.0**(ii-1)) )/(ii-1)

        h += a[1]*numpy.log(T/300)
        h += a[0]*(1/300 - 1/T)

        h += 1.4*287.035*300/(1.4-1)

        return h

    def energy(self, T):

        h = self.enthalpy(T)

        e = h - 287.035*T

        return e

    def temperature(self, e):

        T0 = e/(287.035/(1.4-1))
        e0 = self.energy(T0)

        T1 = T0*1.001
        e1 = self.energy(T1)        

        T2 = T1 + (e - e1)*(T1-T0)/(e1-e0)

        T0 = T1
        e0 = e1
        T1 = T2
        e1 = self.energy(T1)        

        T2 = T1 + (e - e1)*(T1-T0)/(e1-e0)

        T0 = T1
        e0 = e1
        T1 = T2
        e1 = self.energy(T1)        

        T2 = T1 + (e - e1)*(T1-T0)/(e1-e0)

        return T2


if __name__=="__main__":

    g = gas()

    g.calcCoeffAir()

    """

    for ii in range(0, 17):
        T  =300 + ii*100
        print(T, g.cp(T), g.cp2(T))

    for ii in range(0, 17):
        T  =300 + ii*100
        print(T, g.cp(T), g.cv(T), g.cp(T)/g.cv(T))

    for ii in range(0, 17):
        T  =300 + ii*100
        print(T, g.enthalpy(T), 1.4*287.035*T/(1.4-1))

    """

    for ii in range(0, 17):
        T  =300 + ii*100
        e = g.energy(T)
        print(T, e, g.temperature(e))

