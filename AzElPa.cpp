//
// Created by itsfr on 28/04/2024.
//

#include "AzElPa.h"
#include <cmath>
#include <cstdio>
#include "SAT_Const.h"

Matrix AzElPa(const Matrix& s) {

    double rho = sqrt(s(1, 1) * s(1, 1) + s(1, 2) * s(1, 2));

    double Az = atan2(s(1,1), s(1,2));

    if (Az < 0.0) {
        Az += Constants::pi2;
    }

    double El = atan(s(1,3) / rho);

    Matrix dAds(1, 3);
    dAds(1, 1) = s(1, 2) / (rho * rho);
    dAds(1, 2) = -s(1, 1) / (rho * rho);
    dAds(1, 3) = 0.0;

    Matrix dEds(1, 3);
    dEds(1, 1) = -s(1, 1) * s(1, 3) / rho;
    dEds(1, 2) = -s(1, 2) * s(1, 3) / rho;
    dEds(1, 3) = rho;
    double s_dot_s = s(1,1) * s(1,1) + s(1,2) * s(1,2) + s(1,3) * s(1,3);
    dEds = dEds / s_dot_s;


    Matrix result(1, 8);
    result(1,1) = Az;
    result(1,2) = El;
    result(1,3) = dAds(1, 1);
    result(1,4) = dAds(1, 2);
    result(1,5) = dAds(1, 3);
    result(1,6) = dEds(1, 1);
    result(1,7) = dEds(1, 2);
    result(1,8) = dEds(1, 3);
    return result;
}