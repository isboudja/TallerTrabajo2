//
// Created by itsfr on 02/05/2024.
//

#include "math.h"
#include "Mjday_TDB.h"

/**
 * @brief Calcula la Fecha Juliana Modificada en Tiempo Dinámico Bariocéntrico  a partir de la Fecha Juliana Modificada en Tiempo Terrestre  dada.
 *
 * @param Mjd_TT La Fecha Juliana Modificada en Tiempo Terrestre.
 *
 * @return La Fecha Juliana Modificada en Tiempo Dinámico Bariocéntrico.
 */

double Mjday_TDB(double Mjd_TT) {

    double T_TT = (Mjd_TT - 51544.5) / 36525;

    double Mjd_TDB = Mjd_TT + (0.001658 * sin(628.3076 * T_TT + 6.2401) +
                               0.000022 * sin(575.3385 * T_TT + 4.2970) +
                               0.000014 * sin(1256.6152 * T_TT + 6.1969) +
                               0.000005 * sin(606.9777 * T_TT + 4.0212) +
                               0.000005 * sin(52.9691 * T_TT + 0.4444) +
                               0.000002 * sin(21.3299 * T_TT + 5.5431) +
                               0.000010 * sin(628.3076 * T_TT + 4.2490)) / 86400;

    return Mjd_TDB;
}