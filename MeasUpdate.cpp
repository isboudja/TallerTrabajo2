//
// Created by itsfr on 26/05/2024.
//

#include "MeasUpdate.h"

/**
 * @brief el paso de actualización de la medición en un filtro de Kalman.
 *
 * @param x El vector de estado que se actualizará.
 * @param z El valor medido.
 * @param g La función de medición.
 * @param s La desviación estándar del ruido de la medición.
 * @param G La matriz Jacobiana de la función de medición.
 * @param P La matriz de covarianza del estado.
 * @param n El tamaño del vector de estado.
 *
 * @return La matriz de ganancia de Kalman.
 */

Matrix MeasUpdate(Matrix& x, double z,Matrix& g,double s,Matrix& G, Matrix& P, double n) {
    int m = z;
    Matrix Inv_W(m,m);

    for (int i = 0; i < m; ++i) {
        Inv_W(i+1, i+1) = s * s;
    }


    Matrix K = P * G.transpose()/(Inv_W + G * P * G.transpose());


    x = x + K * (g-z);


    Matrix unos(n, n);
    for(int i=1;i<n;i++){
            unos(i,i) = 1.0;
        }

    P = (unos - K * G) * P;

    return K;
}