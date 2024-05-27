//
// Created by itsfr on 02/05/2024.
//

#include "timediff.h"

/**
 * @brief Calcula varias diferencias de tiempo entre diferentes sistemas de tiempo.
 *
 *
 * @param UT1_UTC Diferencia entre UT1 y UTC.
 * @param TAI_UTC Diferencia entre TAI y UTC.
 * @param UT1_TAI Diferencia entre UT1 y TAI (resultado).
 * @param UTC_GPS Diferencia entre UTC y GPS (resultado).
 * @param UT1_GPS Diferencia entre UT1 y GPS (resultado).
 * @param TT_UTC Diferencia entre TT y UTC (resultado).
 * @param GPS_UTC Diferencia entre GPS y UTC (resultado).
 */

void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS, double& UT1_GPS,
              double& TT_UTC, double& GPS_UTC) {
    const double TT_TAI = 32.184;
    const double GPS_TAI = -19.0;
    const double TT_GPS = TT_TAI - GPS_TAI;
    const double TAI_GPS = -GPS_TAI;

    UT1_TAI = UT1_UTC - TAI_UTC;
    const double UTC_TAI = -TAI_UTC;
    UTC_GPS = UTC_TAI - GPS_TAI;
    UT1_GPS = UT1_TAI - GPS_TAI;
    TT_UTC = TT_TAI - UTC_TAI;
    GPS_UTC = GPS_TAI - UTC_TAI;
}