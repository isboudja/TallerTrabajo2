//
// Created by itsfr on 26/05/2024.
//

#include "TimeUpdate.h"

Matrix TimeUpdate( Matrix& P,  Matrix& Phi, double Qdt) {
    return Phi * P * Phi.transpose() + Qdt;
}