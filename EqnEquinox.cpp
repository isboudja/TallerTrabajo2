//
// Created by isboudja on 15/05/2024.
//

#include <valarray>
#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"

/**
 * @brief Calcula la ecuación de los equinoccios para una fecha dada en Tiempo Terrestre (TT) Modificado en formato de Día Juliano (Mjd_TT).
 *
 * @param Mjd_TT Fecha en Tiempo Terrestre (TT) Modificado en formato de Día Juliano Modificado.
 * @return double La ecuación de los equinoccios (radianes).
 */

double EqnEquinox(double Mjd_TT) {
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);
    double mean_obliquity = MeanObliquity(Mjd_TT);
    double EqE = dpsi * cos(mean_obliquity);
    return EqE;
}