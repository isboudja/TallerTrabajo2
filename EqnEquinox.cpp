//
// Created by isboudja on 15/05/2024.
//

#include <valarray>
#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"


double EqnEquinox(double Mjd_TT) {
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);
    double mean_obliquity = MeanObliquity(Mjd_TT);
    double EqE = dpsi * cos(mean_obliquity);
    return EqE;
}