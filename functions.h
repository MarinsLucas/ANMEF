#include<iostream>
#include<cmath>

#define f long double


#define DIRICHLET 0
#define NEUMANN 1

#define EPSLON 0.0001


f cospix(f x)
{
    return cos(M_PI*x);
}

f constante_1(f x)
{
    return 1; 
}

f pi2sinPiX(f x)
{
    return M_PI*M_PI*sin(M_PI*x);
}

f sinpix(f x)
{
    return sin(M_PI * x);
}

f picospix(f x)
{
    return M_PI*cos(M_PI*x);
}

f mpicospix(f x)
{
    return - M_PI*cos(M_PI*x);
}

f pi2cospix(f x)
{
    return M_PI*M_PI*cos(M_PI*x);
}

f solexata2q(f x)
{
    f epslon = EPSLON; 
    f c2 = (pow(M_E, -1.0/sqrt(epslon)) - 1.0)/(pow(M_E, 1.0/sqrt(epslon)) - pow(M_E, -1.0/sqrt(epslon)));
    f c1 = - 1.0 - c2;
    return c1*pow(M_E, -x/sqrt(epslon)) + c2*pow(M_E, x/sqrt(epslon)) + 1.0;
}

f pisinpix(f x)
{
    return M_PI*sin(M_PI*x);
}