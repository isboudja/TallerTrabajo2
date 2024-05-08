//
// Created by itsfr on 01/05/2024.
//

#include <limits>
#include <stdexcept>
#include "Geodetic.h"
#include "SAT_Const.h"

Matrix Geodetic(const Matrix& r) {


    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double epsRequ = eps * Constants::R_Earth;

    double R_equ = Constants::R_Earth;
    double f = Constants::f_Earth;

    double e2 = f * (2.0 - f);
    double X = r(1,1);
    double Y = r(1,2);
    double Z = r(1,3);
    double rho2 = X * X + Y * Y;

    if (std::sqrt(rho2 + Z * Z) == 0.0) {
        throw std::runtime_error("Invalid input in Geodetic constructor");
    }

    double dZ = e2 * Z;
    double Nh;
    double ZdZ;
    double N;
    while (true) {
        ZdZ = Z + dZ;Nh = sqrt(rho2 + ZdZ * ZdZ);
        double SinPhi = ZdZ / Nh;
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        double dZ_new = N * e2 * SinPhi;

        if (abs(dZ - dZ_new) < epsRequ) {
            break;
        }

        dZ = dZ_new;
    }

    double lon = std::atan2(Y, X);
    double lat = std::atan2(ZdZ, std::sqrt(rho2));
    double h = Nh - N;
    double res[] = { lon, lat, h };
    Matrix result = Matrix(1,3,res,3);
    return result;
}