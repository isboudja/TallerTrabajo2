//
// Created by isboudja on 09/05/2024.
//

#include <cmath>
#include "gmst.h"

/**
 * @brief Calcula el Tiempo Sidéreo Medio de Greenwich para una Fecha Juliana Modificada  dada.
 *
 * @param Mjd_UT1 Fecha Juliana Modificada.
 * @return double El Tiempo Sidéreo Medio de Greenwich en radianes.
 */

double gmst(double Mjd_UT1) {
    double Secs = 86400.0;
    double MJD_J2000 = 51544.5;

    double Mjd_0 = std::floor(Mjd_UT1);
    double UT1 = Secs * (Mjd_UT1 - Mjd_0);
    double T_0 = (Mjd_0 - MJD_J2000) / 36525.0;
    double T = (Mjd_UT1 - MJD_J2000) / 36525.0;

    double gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1
                  + (0.093104 - 6.2e-6 * T) * T * T;

    return 2 * M_PI * ((gmst / Secs)-std::floor(gmst / Secs));
}