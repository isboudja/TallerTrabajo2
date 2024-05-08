//
// Created by itsfr on 01/05/2024.
//

#ifndef UNTITLED_SAT_CONST_H
#define UNTITLED_SAT_CONST_H
#include <cmath>

struct Constants {
    // Mathematical constants
    static const double pi2;
    static const double Rad;
    static const double Deg;
    static const double Arcs;

    // General
    static const double MJD_J2000;
    static const double T_B1950;
    static const double c_light;
    static const double AU;

    // Physical parameters of the Earth, Sun and Moon
    static constexpr double R_Earth = 6378.1363e3;
    static const double f_Earth;
    static const double R_Sun;
    static const double R_Moon;

    // Earth rotation
    static const double omega_Earth;

    // Gravitational coefficients
    static const double GM_Earth;
    static const double GM_Sun;
    static const double GM_Moon;
    static const double GM_Mercury;
    static const double GM_Venus;
    static const double GM_Mars;
    static const double GM_Jupiter;
    static const double GM_Saturn;
    static const double GM_Uranus;
    static const double GM_Neptune;
    static const double GM_Pluto;
};



#endif //UNTITLED_SAT_CONST_H
