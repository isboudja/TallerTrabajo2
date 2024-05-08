//
// Created by itsfr on 28/04/2024.
//

#include "Matrix.h"
#include <math.h>


Matrix R_z(double alpha) {

    double C = cos(alpha);
    double S = sin(alpha);

    Matrix rotmat(3,3);

    rotmat(1,1) = C;  rotmat(1,2) =    S;  rotmat(1,3) = 0.0;
    rotmat(2,1) = -1.0*S;  rotmat(2,2) =      C;  rotmat(2,3) =  0.0;
    rotmat(3,1) = 0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) =  1.0;
    return rotmat;
}
