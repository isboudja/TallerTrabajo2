//
// Created by isboudja on 15/05/2024.
//

#include "NutMatrix.h"

#include <cmath>
#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

// Función para calcular la matriz de nutación
Matrix NutMatrix(double Mjd_TT) {
    // Calcular la oblicuidad media de la eclíptica
    double eps = MeanObliquity(Mjd_TT);

    // Calcular la nutación en longitud y oblicuidad
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);

    // Transformación de ecuador y equinoccio medio a verdadero
    Matrix NutMat = R_x(-eps - deps) * R_z(-dpsi) * R_x(+eps);

    return NutMat;
}