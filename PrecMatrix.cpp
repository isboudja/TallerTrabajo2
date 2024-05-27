//
// Created by isboudja on 15/05/2024.
//

#include "PrecMatrix.h"
#include "R_z.h"
#include "R_y.h"
#include "SAT_Const.h"



/**
 * @brief Calcula la matriz de precesión para un intervalo de tiempo dado.
 *
 * @param Mjd_1 Momento inicial del intervalo de tiempo en Tiempo Terrestre Modificado (TTMJD).
 * @param Mjd_2 Momento final del intervalo de tiempo en Tiempo Terrestre Modificado (TTMJD).
 * @return Matriz de precesión que transforma las coordenadas de un sistema de referencia
 *         a otro sistema de referencia debido a la precesión de los ejes de la Tierra.
 *
 */

Matrix PrecMatrix(double Mjd_1, double Mjd_2) {
    Constants c;
    double T  = (Mjd_1 - c.MJD_J2000) / 36525;
    double dT = (Mjd_2 - Mjd_1) / 36525;


    double zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
    ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / (c.Arcs);
    double z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / (c.Arcs);
    double theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
    ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / (c.Arcs);


    Matrix PrecMat = R_z(-z) * R_y(theta) * R_z(-zeta);

    return PrecMat;
}