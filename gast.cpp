//
// Created by isboudja on 15/05/2024.
//

#include <cmath>
#include "gast.h"
#include "gmst.h"


double gast(double Mjd_UT1) {
    double gmstime = gmst(Mjd_UT1); // Calcular GMST
    double eqnequinox = EqnEquinox(Mjd_UT1); // Calcular la ecuación de equinoccio
    double gstime = fmod(gmstime + eqnequinox, 2 * M_PI); // Calcular GAST
    return gstime;
}