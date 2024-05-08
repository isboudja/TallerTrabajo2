//
// Created by itsfr on 02/05/2024.
//

#include "unit.h"
#include "math.h"
Matrix unit(Matrix& vec) {
    const double small = 0.000001;
    double magv = 0.0;
    Matrix outvec(1,3);

    // Calculate the magnitude of the vector
    for (int i = 0; i < 3; i++) {
        magv += vec(1,i+1) * vec(1,i+1);
    }
    magv = sqrt(magv);

    if (magv > small) {
        for (int i = 0; i < 3; i++) {
            outvec(1,i+1) = vec(1,i+1) / magv;
        }
    } else {
        for (int i = 0; i < 3; i++) {
            outvec(1,i+1) = 0.0;
        }
    }
    return outvec;
}