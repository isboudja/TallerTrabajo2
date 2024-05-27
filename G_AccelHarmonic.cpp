//
// Created by itsfr on 25/05/2024.
//

#include "G_AccelHarmonic.h"
#include "AccelHarmonic.h"

/**
 * @brief Calcula el gradiente de la aceleración debido a fuerzas gravitacionales armónicas.
 *
 * @param r Vector de posición (Matriz 3x1) en el marco inercial.
 * @param U Vector de coeficientes del campo de gravedad.
 * @param n_max Grado máximo del campo de gravedad.
 * @param m_max Orden máximo del campo de gravedad.
 *
 * @return Matrix La matriz de gradiente (Matriz 3x3) que representa las derivadas parciales de los componentes de la aceleración
 */

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    double d = 1.0;
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

            for(l=0;l<3;l++){
                dr(l+1,1) = 0.0;
            }
                dr(i+1,1) = d;

         r_plus_dr= r;
         r_minus_dr = r;
        for (int j = 0; j < 3; ++j) {
            r_plus_dr(j+1,1) += dr(j+1,1) / 2;
            r_minus_dr(j+1,1) -= dr(j+1,1) / 2;
        }



        accel_plus = AccelHarmonic(r_plus_dr, U, n_max, m_max);
        accel_minus = AccelHarmonic(r_minus_dr, U, n_max, m_max);

        for (int j = 0; j < 3; ++j) {
            da(j+1,1) = (accel_plus(1,j+1) - accel_minus(1,j+1)) / d;
        }


        for (int j = 0; j < 3; ++j) {
            G(j+1,i+1) = da(j+1,1);
        }
    }
    return G;
}