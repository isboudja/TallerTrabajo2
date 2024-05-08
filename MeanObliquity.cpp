//
// Created by itsfr on 02/05/2024.
//

#include "MeanObliquity.h"
#include "SAT_Const.h"
#include <cmath>

const double MJD_J2000 = 51544.5; // Modifica según el valor de const.MJD_J2000 en tu código
const double Rad = M_PI / 180.0;

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - Constants::MJD_J2000) / 36525.0;

    double MOblq = Constants::Rad * (84381.448 / 3600 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600);

    return MOblq;
}