//
// Created by itsfr on 28/04/2024.
//

#ifndef UNTITLED_CHEB3D_H
#define UNTITLED_CHEB3D_H


#include <cmath>
#include <stdexcept>
#include "Matrix.h"


Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);


#endif //UNTITLED_CHEB3D_H
