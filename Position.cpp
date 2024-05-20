//
// Created by itsfr on 02/05/2024.
//

#include "Position.h"
#include "SAT_Const.h"
void Position(double lon, double lat, double h, Matrix& r) {

    double R_equ = Constants::R_Earth;
    double f = Constants::f_Earth;

    double e2 = f * (2.0 - f);
    double CosLat = cos(lat);
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    r(1,1) = (N + h) * CosLat * cos(lon);
    r(1,2) = (N + h) * CosLat * sin(lon);
    r(1,3) = ((1.0 - e2) * N + h) * SinLat;
}