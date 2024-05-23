//
// Created by itsfr on 01/05/2024.
//

#include "Legendre.h"

#include <cmath>
#include "Matrix.h" // Asegúrate de incluir la definición de la clase Matrix

Matrix Legendre(int n, int m, double fi) {

    Matrix pnm(1, (n + 1)*(m + 1));
    Matrix dpnm(1, (n + 1)*(m + 1));
    Matrix result(1,(n+1)*(m+1)*2);
    pnm(1,1)=1.0;
    dpnm(1,1)=0.0;
    pnm(1,5)=sqrt(3.)*cos(fi);
    dpnm(1,5)=-sqrt(3.)*sin(fi);

    int g = (n+1)*2;
    g = g - (n-1);
    for (int i = 2; i <= n; i++) {
        pnm(1, g+n+2) = sqrt((2. * i + 1.) / (2.0 * i)) * cos(fi) * pnm(1, g);
        g= g+(n+2);
    }

    g = (n+1)*2;
    g = g -(n-1);
    for (int i = 2; i <= n; i++) {
        dpnm(1, g+n+2) = sqrt((2. * i + 1.) / (2.0 * i)) * ((cos(fi) * dpnm(1, g)) - (sin(fi) * pnm(1, g)));
        g= g+(n+2);
    }

    g = (n+2);
    int h = g+n;
    for (int i = 1; i <= n; i++) {
        pnm(1, g) = sqrt(2. * i + 1.) * sin(fi) * pnm(1,h-g-1);
        g= g+(n+2);
        h = h+g;
    }

    g = (n+2);
    h = g +n;
    for (int i = 1; i <= n; i++) {
        dpnm(1, g) = sqrt(2. * i + 1.)*((cos(fi) * pnm(1, h-g-1))+(sin(fi) * dpnm(1, h-g-1)));
        g= g+(n+2);
        h = h+g;
    }

    double j = 0., k = 2.;
    g = (n+2);
    h = g+1+n;
    while (true) {
        for (int i = k; i <= n; ++i) {
            pnm(1,h) = sqrt((2. * i + 1.) / ((i - j) * (i + j))) *
                       ((sqrt(2 * i - 1) * sin(fi) * pnm(1, g)) -
                        (sqrt(((i + j - 1.) * (i - j - 1.)) / (2. * i - 3.)) * pnm(1, h-(n+1)*2)));

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
            dpnm(1,h)= sqrt((2.*i+1.)/((i-j)*(i+j)))*((sqrt(2.*i-1.)*sin(fi)*dpnm(1.,g))
                    +(sqrt(2.*i-1.)*cos(fi)*pnm(1,g))-(sqrt(((i+j-1.)*(i-j-1.))/(2.*i-3.))*dpnm(1,h-(n+1)*2)));


        }
        ++j;
        ++k;
        h = h+n+1;
        g = g+n+1;
        if (j > m) {
            break;
        }
    }

    for(h = 0;h<(n+1)*(m+1);h++){
        result(1,h+1) = pnm(1, h + 1);
    }
    g = 0;
    for(h = (n+1)*(m+1);h<(n+1)*(m+1)*2;h++){
        result(1,h+1) = dpnm(1, g + 1);
        g++;
    }
    return result;

}