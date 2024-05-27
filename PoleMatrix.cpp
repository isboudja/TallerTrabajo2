//
// Created by isboudja on 15/05/2024.
//

#include "PoleMatrix.h"
#include "R_y.h"
#include "R_x.h"

/**
 * @brief Calcula la matriz de rotación asociada a las coordenadas de los polos celestes para una fecha dada.
 *
 * @param xp La coordenada polar celeste en el eje y en radianes.
 * @param yp La coordenada polar celeste en el eje x en radianes.
 * @return La matriz de rotación asociada a las coordenadas de los polos celestes.
 *
 */

Matrix PoleMatrix(double xp, double yp) {
    Matrix PoleMat = R_y(-xp) * R_x(-yp);

    return PoleMat;
}