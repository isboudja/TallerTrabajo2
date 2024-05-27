//
// Created by itsfr on 01/05/2024.
//

#include <cmath>
#include "Frac.h"

/**
 * @brief Calcula la parte fraccionaria de un número de precisión doble.
 *
 * @param x El número de precisión doble de entrada.
 * @return double La parte fraccionaria del número de entrada.
 */

double Frac(double x) {
    return x - std::floor(x);
}