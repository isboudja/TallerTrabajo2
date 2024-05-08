//
// Created by itsfr on 02/05/2024.
//

#include "sign_.h"
#include "math.h"
double sign_(double a, double b) {
    if (b >= 0.0) {
        return abs(a);
    } else {
        return -abs(a);
    }
}