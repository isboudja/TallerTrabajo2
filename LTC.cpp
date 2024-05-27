//
// Created by itsfr on 26/05/2024.
//

#include "LTC.h"
#include "R_y.h"
#include "R_z.h"

/**
 * Calcula la matriz de transformación de coordenadas geodésicas (longitud, latitud)
 * a coordenadas cartesianas (x, y, z) utilizando el método de Coordenadas Tangentes Locales (LTC).
 *
 * @param lon La longitud en radianes.
 * @param lat La latitud en radianes.
 *
 * @return La matriz de transformación de coordenadas geodésicas a coordenadas cartesianas.
 */

Matrix LTC(double lon, double lat) {
    Matrix M = R_y(-1.0 * lat) * R_z(lon);

    for (int j = 0; j < 3; ++j) {
        double aux = M(1, j+1);
        M(1, j+1) = M(2, j+1);
        M(2, j+1) = M(3, j+1);
        M(3, j+1) = aux;
    }

    return M;
}