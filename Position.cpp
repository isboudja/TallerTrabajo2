//
// Created by itsfr on 02/05/2024.
//

#include "Position.h"
#include "SAT_Const.h"

/**
 * @brief Calcula la posición tridimensional (coordenadas cartesianas) de un punto sobre la Tierra.
 *
 * @param lon Longitud del punto en radianes.
 * @param lat Latitud del punto en radianes.
 * @param h Altura sobre el elipsoide de referencia en metros.
 * @param r Matriz de 1x3 donde se almacenarán las coordenadas cartesianas resultantes [x, y, z].
 *
 */
void Position(double lon, double lat, double h, Matrix& r) {

    double R_equ = Constants::R_Earth;
    double f = Constants::f_Earth;

    double e2 = f * (2.0 - f);
    double CosLat = cos(lat);
    double SinLat = sin(lat);


    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    r(1,1) = (N + h) * CosLat * cos(lon);
    r(1,2) = (N + h) * CosLat * sin(lon);
    r(1,3) = ((1.0 - e2) * N + h) * SinLat;
}