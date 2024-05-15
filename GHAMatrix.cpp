//
// Created by isboudja on 15/05/2024.
//

#include "GHAMatrix.h"
#include "gast.h"
#include "R_z.h"

Matrix GHAMatrix(double Mjd_UT1) {

    double GAST = gast(Mjd_UT1);


    Matrix GHAmat = R_z(GAST);

    return GHAmat;
}