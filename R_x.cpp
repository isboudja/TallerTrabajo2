


#include "Matrix.h"
#include <math.h>

/**
 * @brief la matriz de rotación alrededor del eje X.
 *
 * @param alpha Ángulo de rotación en radianes.
 * @return Matriz de rotación 3x3 que rota los vectores alrededor del eje X
 *         por el ángulo especificado.
 *
 */

Matrix R_x(double alpha) {

    double C = cos(alpha);
    double S = sin(alpha);

    Matrix rotmat(3,3);

    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;
    return rotmat;
}