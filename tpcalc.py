import numpy
import matplotlib.pyplot as plt


class gas():

    def __init__(self):

        """
        Based on:

        Witte, D. W. and Tatum K. E., Computer Code for Determination Perfect Gas Properties of Thermally Perfect Gas Properties, 1994


        """
        
        self.R = 8.31446261815324
        self.fracN2 = 0.75530
        self.fracO2 = 0.23140
        self.fracAr = 0.01290
        self.fracCO2 = 0.00040
        self.Rgas = 287.035

    def cpN2(self, T):

        R = self.R/0.0280160

        if T < 1000:
            a = [-1.33984200E+01, 1.34280300E+00, 3.45742000E+00, 5.74727600E-04,-3.21711900E-06, 7.50775400E-09, -5.90150500E-12, 1.50979900E-15]

        else:
            a = [5.87702841E+05, -2.23921563E+03, 6.06686971E+00, -6.13957913E-04, 1.49178026E-07, -1.92307130E-11, 1.06193594E-15, 0.00000000E+00]

        cp = 0.0        
        for ii in range(0, 8):
            cp += a[ii]*(T**(ii-2))

        cp *= R

        return cp

    def cpO2(self, T):

        R = self.R/0.032

        if T < 1000:
            a = [3.88517500E+01, -2.70630800E+00, 3.56119600E+00, -3.32782400E-04, -1.18148000E-06, 1.10853500E-08, -1.49299400E-11, 5.99553800E-15]

        else:
            a = [-1.05642070E+06, 2.41123849E+03, 1.73474238E+00, 1.31512292E-03,-2.29995151E-07, 2.13144378E-11, -7.87498771E-16, 0.00000000E+00]

        cp = 0.0        
        for ii in range(0, 8):
            cp += a[ii]*(T**(ii-2))

        cp *= R

        return cp

    def cpAr(self, T):

        R = self.R/0.039944
        cp = 2.5*R

        return cp

    def cpCO2(self, T):

        R = self.R/0.040220

        if T < 1000:
            a = [-5.34531900E+03, 2.28785400E+02, 2.56174500E-01, 1.68213800E-02,-2.09949800E-05, 1.40430800E-08, -3.81732100E-12, 0.00000000E+00]

        else:
            a = [1.15460081E+05, -1.78337370E+03, 8.28644239E+00,  -8.98356945E-05, 4.26107946E-09,  -1.81443266E-12,  6.29130739E-16, 0.00000000E+00]

        cp = 0.0        
        for ii in range(0, 8):
            cp += a[ii]*(T**(ii-2))

        cp *= R

        return cp


    def cp(self, T):

        return self.fracN2*self.cpN2(T) + self.fracO2*self.cpO2(T) + self.fracAr*self.cpAr(T) + self.fracCO2*self.cpCO2(T)

    def cv(self, T):

        return self.cp(T) - self.Rgas

if __name__=="__main__":

    g = gas()

    for ii in range(0, 17):
        T  =300 + ii*100
        print(T, g.cp(T), g.cv(T), g.cp(T)/g.cv(T))
