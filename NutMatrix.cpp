//
// Created by isboudja on 15/05/2024.
//

#include "NutMatrix.h"

#include <cmath>
#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

/**
 * @brief Calcula la matriz de rotación asociada a las variaciones en la oblicuidad eclíptica y en la longitud eclíptica para una fecha dada en Tiempo Terrestre (TT).
 *
 * @param Mjd_TT La Fecha Juliana Modificada (MJD) en Tiempo Terrestre (TT).
 * @return La matriz de rotación asociada a las variaciones en la oblicuidad eclíptica y en la longitud eclíptica.
 *
 */

Matrix NutMatrix(double Mjd_TT) {

    double eps = MeanObliquity(Mjd_TT);
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);
    Matrix NutMat = R_x(-eps - deps) * R_z(-dpsi) * R_x(+eps);

    return NutMat;
}