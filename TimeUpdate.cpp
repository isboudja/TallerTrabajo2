//
// Created by itsfr on 26/05/2024.
//

#include "TimeUpdate.h"

/**
 *
 * @param P La matriz de covarianza del error de estado en el tiempo anterior.
 * @param Phi La matriz de propagación del estado.
 * @param Qdt La matriz de covarianza del proceso escalada por el tiempo.
 * @return La matriz de covarianza del error de estado actualizada en el tiempo.
 */

Matrix TimeUpdate( Matrix& P,  Matrix& Phi, double Qdt) {
    return Phi * P * Phi.transpose() + Qdt;
}