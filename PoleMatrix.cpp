//
// Created by isboudja on 15/05/2024.
//

#include "PoleMatrix.h"
#include "R_y.h"
#include "R_x.h"

Matrix PoleMatrix(double xp, double yp) {
    Matrix PoleMat = R_y(-xp) * R_x(-yp);

    return PoleMat;
}