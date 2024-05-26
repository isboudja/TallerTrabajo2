//
// Created by itsfr on 25/05/2024.
//

#include "G_AccelHarmonic.h"
#include "AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max,Matrix &Snm,Matrix &Cnm) {
    double d = 1.0;   // Position increment [m]
    Matrix G(3, 3);
    Matrix dr(3, 1.0);
    Matrix r_plus_dr(3,1);
    Matrix r_minus_dr(3,1);
    Matrix accel_plus(1,3);
    Matrix da(3, 1.0);
    Matrix accel_minus(1,3);
    // Gradient
    int i;
    int l;
    for (i = 0; i < 3; i++) {
        // Set offset in i-th component of the position vector
            for(l=0;l<3;l++){
                dr(l+1,1) = 0.0;
            }
                dr(i+1,1) = d;
        // Create r+dr/2 and r-dr/2 vectors
         r_plus_dr= r;
         r_minus_dr = r;
        for (int j = 0; j < 3; ++j) {
            r_plus_dr(j+1,1) += dr(j+1,1) / 2;
            r_minus_dr(j+1,1) -= dr(j+1,1) / 2;
        }

        // Acceleration difference

        accel_plus = AccelHarmonic(r_plus_dr, U, n_max, m_max,Snm,Cnm);
        accel_minus = AccelHarmonic(r_minus_dr, U, n_max, m_max,Snm,Cnm);

        for (int j = 0; j < 3; ++j) {
            da(j+1,1) = (accel_plus(1,j+1) - accel_minus(1,j+1)) / d;
        }

        // Derivative with respect to i-th axis
        for (int j = 0; j < 3; ++j) {
            G(j+1,i+1) = da(j+1,1);
        }
    }
    return G;
}