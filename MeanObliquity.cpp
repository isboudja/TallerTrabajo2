//
// Created by itsfr on 02/05/2024.
//

#include "MeanObliquity.h"
#include "SAT_Const.h"
#include <cmath>

const double MJD_J2000 = 51544.5;
const double Rad = M_PI / 180.0;

/**
 * Calcula la oblicuidad media de la eclíptica para una fecha juliana de Tiempo Terrestre dada.
 *
 * @param Mjd_TT La fecha juliana de Tiempo Terrestre .
 *
 * @return La oblicuidad media de la eclíptica en radianes.
 */

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - Constants::MJD_J2000) / 36525.0;

    double MOblq = Constants::Rad * (84381.448 / 3600 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600);

    return MOblq;
}