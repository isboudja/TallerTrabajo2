//
// Created by isboudja on 15/05/2024.
//

#include <cmath>
#include "gast.h"
#include "gmst.h"
#include "EqnEquinox.h"
/**
 * @brief Calcula el Tiempo Sidéreo Aparente de Greenwich.
 *
 * @param Mjd_UT1 Fecha Juliana Modificada.
 * @return double El Tiempo Sidéreo Aparente de Greenwich en radianes.
 */

double gast(double Mjd_UT1) {
    double gmstime = gmst(Mjd_UT1);
    double eqnequinox = EqnEquinox(Mjd_UT1);
    double gstime = fmod(gmstime + eqnequinox, 2 * M_PI);
    return gstime;
}