//
// Created by itsfr on 26/05/2024.
//

#include "MeasUpdate.h"

Matrix MeasUpdate(Matrix& x, double z,Matrix& g,double s,Matrix& G, Matrix& P, double n) {
    int m = z;
    Matrix Inv_W(m,m);
    // Calculate Inverse weight (measurement covariance)
    for (int i = 0; i < m; ++i) {
        Inv_W(i+1, i+1) = s * s;
    }

    // Kalman gain
    Matrix K = P * G.transpose()/(Inv_W + G * P * G.transpose());

    // State update
    x = x + K * (g-z);

    // Covariance update
    Matrix unos(n, n);
    for(int i=1;i<n;i++){
            unos(i,i) = 1.0;
        }

    P = (unos - K * G) * P;

    return K;
}