//
// Created by itsfr on 02/05/2024.
//

#include <cmath>
#include "Mjday.h"

/**
 * @brief Calcula la Fecha Juliana Modificada  a partir de los componentes de fecha y tiempo dados.
 *
 * @param yr  El año.
 * @param mon El mes.
 * @param day El día del mes.
 * @param hr  La hora del día.
 * @param min El minuto de la hora.
 * @param sec El segundo del minuto.
 *
 * @return La Fecha Juliana Modificada.
 */

double Mjday(double yr, double mon, double day, double hr, double min, double sec) {
    double jd = 367.0 * yr - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25) +
                floor(275 * mon / 9.0) +
                day + 1721013.5 +
                ((sec / 60.0 + min) / 60.0 + hr) / 24.0;

    double Mjd = jd - 2400000.5;

    return Mjd;
}