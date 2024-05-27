//
// Created by itsfr on 01/05/2024.
//

#include "Legendre.h"

#include <cmath>
#include "Matrix.h" // Asegúrate de incluir la definición de la clase Matrix

/**
 * @brief las funciones Legendre asociadas y sus derivadas hasta el grado n y orden m,evaluadas en la latitud dada.
 *
 * @param n El grado máximo de los polinomios de Legendre.
 * @param m El orden máximo de los polinomios de Legendre.
 * @param fi La latitud en radianes donde se evalúan las funciones Legendre asociadas.
 *
 * @return Una matriz que contiene las funciones Legendre asociadas y sus derivadas
 */

Matrix Legendre(int n, int m, double fi) {

    Matrix pnm((n + 1), (m + 1));
    Matrix dpnm((n + 1), (m + 1));
    Matrix result(1,(n+1)*(m+1)*2);
    pnm(1,1)=1.0;
    dpnm(1,1)=0.0;
    pnm(2,2)=sqrt(3.)*cos(fi);
    dpnm(2,2)=-sqrt(3.)*sin(fi);

    int g = (n + 1)*(m + 1);
    for (int i = 2; i <= n; i++) {
        pnm(i+1, i+1) = sqrt((2. * i + 1.) / (2.0 * i)) * cos(fi) * pnm(i, i);
        g= g+(n+1);
    }

    g = n+2;
    for (int i = 2; i <= n; i++) {
        dpnm(i+1, i+1) = sqrt((2. * i + 1.) / (2.0 * i)) * ((cos(fi) * dpnm(i, i)) - (sin(fi) * pnm(i, i)));
        g= g+(n+1);
    }

        g = n+2;
    for (int i = 1; i <= n; i++) {
        pnm(i+1, i) = sqrt(2. * i + 1.) * sin(fi) * pnm(i,i);
        g= g+(n+2);
    }

    g = (n+2);
    for (int i = 1; i <= n; i++) {
        dpnm(i+1, i) = sqrt(2. * i + 1.)*((cos(fi) * pnm(i, i))+(sin(fi) * dpnm(i, i)));
        g= g+n+2;
    }

    double j = 0., k = 2.;
    g = (n+2);
    int h = g+1+n;
    while (true) {
        for (int i = k; i <= n; ++i) {
            pnm(i+1,j+1) = sqrt((2. * i + 1.) / ((i - j) * (i + j))) *
                       ((sqrt(2 * i - 1) * sin(fi) * pnm(i,j+1)) -
                        (sqrt(((i + j - 1.) * (i - j - 1.)) / (2. * i - 3.)) * pnm(i-1,j+1)));


        }
        ++j;
        ++k;
        h = h+n+1;
        g = g+n+1;
        if (j > m) {
            break;
        }
    }
    j = 0., k = 2.;
    g = (n+2);
    h = g+1+n;
    while (true) {
        for (int i = k; i <= n; ++i) {
            dpnm(i+1,j+1)= sqrt((2.*i+1.)/((i-j)*(i+j)))*((sqrt(2.*i-1.)*sin(fi)*dpnm(i,j+1))
                    +(sqrt(2.*i-1.)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1.)*(i-j-1.))/(2.*i-3.))*dpnm(i-1,j+1)));



        }
        ++j;
        ++k;
        h = h+n+1;
        g = g+n+1;
        if (j > m) {
            break;
        }
    }
    g = 1;
    for(j=0;j<pnm.fil;j++){
        for(h = 0;h<pnm.col;h++){
            result(1,g) = pnm(j+1, h + 1);
            g++;
        }
    }

    for(j=0;j<dpnm.fil;j++){
        for(h = 0;h<dpnm.col;h++){
            result(1,g) = dpnm(j+1, h + 1);
            g++;
        }
    }

    return result;

}