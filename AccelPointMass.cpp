//
// Created by itsfr on 28/04/2024.
//

#include <cstdio>
#include "AccelPointMass.h"
#include "math.h"

/**
 * @brief Calcula la aceleración debido a una masa puntual.
 *
 * @param r Una referencia a una Matriz que representa el vector de posición del satélite (matriz 3x1).
 * @param s Una referencia a una Matriz que representa el vector de posición de la masa puntual (matriz 3x1).
 * @param GM El parámetro gravitacional (parámetro gravitacional estándar) de la masa puntual.
 *
 * @return Una Matriz que contiene el vector de aceleración (matriz 3x1).
 */

Matrix AccelPointMass(const Matrix& r, const Matrix& s, double GM) {

    Matrix d = r-s;


    double norm_d = d.norm();
    double norm_s = s.norm();


    return (d/pow(norm_d,3) + s/pow(norm_s,3))*-GM;
}