//
// Created by itsfr on 28/04/2024.
//

#include "AzElPa.h"
#include <cmath>
#include <cstdio>
#include "SAT_Const.h"

/**
 * @brief Calcula el azimut, la elevación y sus derivadas parciales.
 *
 * @param s Una referencia a una Matriz que representa el vector de posición (matriz 1x3).
 *
 * @return Una Matriz que contiene el azimut, la elevación y sus derivadas parciales (matriz 1x8):
 * - result(1, 1): Ángulo de azimut (Az) en radianes.
 * - result(1, 2): Ángulo de elevación (El) en radianes.
 * - result(1, 3): Derivada parcial de Az respecto a x.
 * - result(1, 4): Derivada parcial de Az respecto a y.
 * - result(1, 5): Derivada parcial de Az respecto a z.
 * - result(1, 6): Derivada parcial de El respecto a x.
 * - result(1, 7): Derivada parcial de El respecto a y.
 * - result(1, 8): Derivada parcial de El respecto a z.
 */

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