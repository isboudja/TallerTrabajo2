//
// Created by isboudja on 15/05/2024.
//

#include "GHAMatrix.h"
#include "gast.h"
#include "R_z.h"

/**
 * @brief Calcula la matriz del Ángulo Horario de Greenwich.
 *
 *
 * @param Mjd_UT1 Fecha Juliana Modificada.
 * @return Matrix La matriz del Ángulo Horario de Greenwich.
 */

Matrix GHAMatrix(double Mjd_UT1) {

    double GAST = gast(Mjd_UT1);


    Matrix GHAmat = R_z(GAST);

    return GHAmat;
}