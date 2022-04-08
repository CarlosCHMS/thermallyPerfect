#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct
{
    
    int type;
    double Rgas;
    double* a0;
    double* a1;    
            
} GAS;

GAS* gasInit()
{

    GAS* gas = malloc(sizeof(GAS));
    gas->Rgas = 287.035;

    double a0[8] = {-1.11245334e+03,  1.57330387e+02,  9.95843854e+02,  1.10220174e-01, -7.93915153e-04,  2.35056151e-06, -2.22081406e-09,  6.98903775e-13};
    double a1[8] = { 6.82296802e+07, -3.57105636e+05,  1.47161884e+03, -5.85585389e-02, 1.96110679e-05, -3.02929435e-09,  1.90742602e-13,  0.00000000e+00};

    gas->a0 = malloc(8*sizeof(double));
    gas->a1 = malloc(8*sizeof(double));
    
    for(int ii=0; ii<8; ii++)
    {
        gas->a0[ii] = a0[ii];
        gas->a1[ii] = a1[ii];        
    }
    
    return gas;

}

double gasEnthalpy(GAS* gas, double T)
{
    double h = 0.0;
    double* a;
    double Ta0;
    double Ta1;

    if(gas->type == 0)
    {
        h = 1.4*gas->Rgas*T/(1.4-1);
    }
    else
    {    

        if(T<1000)
        {
            a = gas->a0;
        }
        else
        {
            a = gas->a1;
        }
                
        Ta0 = 300;
        Ta1 = T;        
        for(int ii=2; ii<8; ii++)
        {
            h += a[ii]*( Ta1 - Ta0 )/(ii-1);
            Ta0 *= 300;
            Ta1 *= T;            
        }

        h += a[1]*log(T/300);
        h += a[0]*(1./300. - 1./T);
        h += 1.4*gas->Rgas*300/(1.4-1.);        
    }

    return h;

}

double gasEnergy(GAS* gas, double T)
{

    return gasEnthalpy(gas, T) - gas->Rgas*T;

}

double gasTemperature(GAS* gas, double e)
{

    double T;
    double T0, T1, T2, e0, e1;

    if(gas->type == 0)
    {
        T = e/(gas->Rgas/(1.4-1.));
    }
    else
    {   
        T0 = e/(gas->Rgas/(1.4-1));
        e0 = gasEnergy(gas, T0);

        T1 = T0*1.1;
        e1 = gasEnergy(gas, T1);

        T2 = T1 + (e - e1)*(T1-T0)/(e1-e0);

        while(fabs(T2 - T1)>1.e-6)
        {
            T0 = T1;
            e0 = e1;
            T1 = T2;
            e1 = gasEnergy(gas, T1);

            T2 = T1 + (e - e1)*(T1-T0)/(e1-e0);
        }
        
        T = T2;
    }

    return T;

}

double gasCp(GAS* gas, double T)
{

    double cp;
    double* a;
    double Ta0;

    if(gas->type == 0)
    {
        cp = 1.4*gas->Rgas/(1.4-1.);
    }
    else
    {    
        if(T<1000)
        {
            a = gas->a0;
        }
        else
        {
            a = gas->a1;
        }
        
        cp = 0.0;
        Ta0 = 1/(T*T);
        for(int ii=0; ii<8; ii++)
        {
            cp += a[ii]*Ta0;
            Ta0 *= T;
        }  
    }

    return cp;

}

double gasSound(GAS* gas, double T)
{

    double cp = gasCp(gas, T);
    double g = cp/(cp - gas->Rgas);
    
    return sqrt(g*gas->Rgas*T);

}

int main()
{

    GAS* gas = gasInit();
    gas->type = 0;

    double T0 = 2000;
    

    double e = gasEnergy(gas, T0);
    double T = gasTemperature(gas, e);
    double Cp = gasCp(gas, T0);
    double a = gasSound(gas, T0);

    printf("%f, %f, %f, %f, %f\n", T0, e, T, Cp, a);

    gas->type = 1;

    e = gasEnergy(gas, T0);
    T = gasTemperature(gas, e);
    Cp = gasCp(gas, T0);
    a = gasSound(gas, T0);

    printf("%f, %f, %f, %f, %f\n", T0, e, T, Cp, a);

    double Cp1 = (gasEnthalpy(gas, T0*1.001) - gasEnthalpy(gas, T0))/(T0*0.001);

    printf("%f, %f\n", Cp, Cp1);

    return 0;

}

