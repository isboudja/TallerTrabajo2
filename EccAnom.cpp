//
// Created by itsfr on 02/05/2024.
//

#include "EccAnom.h"

#include <cmath>
#include <stdexcept>
#include <limits>

double EccAnom(double M, double e) {
    const int m = 15;
    int i = 1;

    M = fmod(M, 2.0 * M_PI);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    while (fabs(f) > 1e2 * std::numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        ++i;
        if (i == m) {
            throw std::runtime_error("problemas en EccAnom");
        }
    }

    return E;
}