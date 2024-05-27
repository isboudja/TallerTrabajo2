//
// Created by itsfr on 28/04/2024.
//

#include "Cheb3D.h"

#include <cmath>
#include <stdexcept>


/**
 * @brief Calcula los valores interpolados utilizando polinomios de Chebyshev en 3 dimensiones.
 *
 * @param t El variable de tiempo para el cual se realiza la interpolación.
 * @param N El orden de los polinomios de Chebyshev.
 * @param Ta El tiempo de inicio del intervalo de interpolación.
 * @param Tb El tiempo de finalización del intervalo de interpolación.
 * @param Cx Los coeficientes de los polinomios de Chebyshev para la dimensión x.
 * @param Cy Los coeficientes de los polinomios de Chebyshev para la dimensión y.
 * @param Cz Los coeficientes de los polinomios de Chebyshev para la dimensión z.
 *
 * @return Una matriz que contiene los valores interpolados en tres dimensiones (matriz 1x3).
 */

Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz){
    if (t < Ta || t > Tb) {
        throw std::runtime_error("ERROR: Tiempo fuera del rango en Cheb3D::Value\n");
    }

    double tau = (2 * t - Ta - Tb) / (Tb - Ta);
    Matrix f1(1, 3);
    Matrix f2(1, 3);
    Matrix old_f1(1,3);
    for (int i = N; i >= 2; --i) {
        old_f1 = f1;
        f1(1,1) = 2*tau*f1(1,1) - f2(1,1) + Cx(1,i);
        f1(1,2) = 2 * tau * f1(1,2) - f2(1,2) + Cy(1,i);
        f1(1,3) = 2 * tau * f1(1,3) - f2(1,3) + Cz(1,i);
        f2 = old_f1;
    }

    Matrix ChebApp(1,3);
    ChebApp(1,1) = tau * f1(1,1)- f2(1,1) + Cx(1,1);
    ChebApp(1,2) = tau * f1(1,2)- f2(1,2) + Cy(1,1);
    ChebApp(1,3) = tau * f1(1,3) - f2(1,3) + Cz(1,1);

    return ChebApp;
}