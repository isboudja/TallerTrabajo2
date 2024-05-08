//
// Created by itsfr on 02/05/2024.
//

#include "timediff.h"

void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS, double& UT1_GPS,
              double& TT_UTC, double& GPS_UTC) {
    const double TT_TAI = 32.184;           // TT-TAI time difference [s]
    const double GPS_TAI = -19.0;           // GPS-TAI time difference [s]
    const double TT_GPS = TT_TAI - GPS_TAI; // TT-GPS time difference [s]
    const double TAI_GPS = -GPS_TAI;        // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC - TAI_UTC;           // UT1-TAI time difference [s]
    const double UTC_TAI = -TAI_UTC;       // UTC-TAI time difference [s]
    UTC_GPS = UTC_TAI - GPS_TAI;           // UTC_GPS time difference [s]
    UT1_GPS = UT1_TAI - GPS_TAI;           // UT1-GPS time difference [s]
    TT_UTC = TT_TAI - UTC_TAI;             // TT-UTC time difference [s]
    GPS_UTC = GPS_TAI - UTC_TAI;           // GPS-UTC time difference [s]
}