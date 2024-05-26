//
// Created by itsfr on 28/04/2024.
//

#include <cstdio>
#include "AccelPointMass.h"
#include "math.h"

Matrix AccelPointMass(const Matrix& r, const Matrix& s, double GM) {
    // Relative position vector of satellite w.r.t. point mass
    Matrix d = r-s;

    // Calculate the norms
    double norm_d = d.norm();
    double norm_s = s.norm();
    // Acceleration

    return (d/pow(norm_d,3) + s/pow(norm_s,3))*-GM;
}