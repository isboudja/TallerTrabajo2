//
// Created by itsfr on 02/05/2024.
//

#include "sign_.h"
#include "math.h"

/**
 * @brief Calcula el valor absoluto de a con el signo de b.
 *
 * @param a El valor cuyo signo se conservará.
 * @param b El valor que determina el signo del resultado.
 * @return El valor absoluto de a con el mismo signo que b.
 *
 */
double sign_(double a, double b) {
    if (b >= 0.0) {
        return abs(a);
    } else {
        return -abs(a);
    }
}