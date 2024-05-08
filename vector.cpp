#include <cmath>
#include "vector.h"

// Función que calcula la norma de un vector
// In: v vector (double)
// In: n dimensión (int)
// Out: devuelve la norma de v
double norm(double v[], int n)
{
	double suma = 0.0;
	int i;
	
	if(n <= 0)
		throw "Empty vector";

	for(i = 0; i < n; ++i)
		suma += v[i]*v[i];

	return(sqrt(suma));
}

// Función que calcula el producto escalar de dos vectores
// In: v1, v2 vectores (double)
// In: n1, n2 dimensiones (int)
// Out: devuelve el producto escalar de v1 y v2
double dot(double v1[], double v2[], int n1, int n2)
{
        double suma = 0.0;
        int i;

        if(n1 <= 0 || n2 <= 0 || n1 !=n2)
                throw "Different dimensions";

        for(i = 0; i < n1; ++i)
                suma += v1[i]*v2[i];

        return(suma);
}

// Función que calcula el producto vectorial de dos vectores
// In: v1, v2 vectores (double)
// In: nv1, nv2 dimensiones (int)
// Out: devuelve el producto vectorial de v1 y v2
void cross(double v[], int &nv, double v1[], double v2[], int nv1, int nv2)
{
     if((nv1 <= 0) || (nv2 <= 0) | (nv1 != nv2))
        throw "Empty vector or different dimensions";

    v[0] = -v1[2]*v2[1] + v1[1]*v2[2];
    v[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v[2] =  -v1[1]*v2[0] + v1[0]*v2[1];
    nv = nv1;
}

