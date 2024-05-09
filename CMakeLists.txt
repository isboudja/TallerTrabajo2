cmake_minimum_required(VERSION 3.23)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled Matrix.cpp Matrix.h TestMain.cpp R_x.h globals.cpp globals.h R_x.cpp R_y.cpp R_y.h R_z.cpp R_z.h AccelPointMass.cpp AccelPointMass.h AzElPa.cpp AzElPa.h Cheb3D.cpp Cheb3D.h SAT_Const.cpp SAT_Const.h Frac.cpp Frac.h Geodetic.cpp Geodetic.h Legendre.cpp Legendre.h EccAnom.cpp EccAnom.h MeanObliquity.cpp MeanObliquity.h Mjday.cpp Mjday.h Mjday_TDB.cpp Mjday_TDB.h NutAngles.cpp NutAngles.h Position.cpp Position.h sign_.cpp sign_.h timediff.cpp timediff.h unit.cpp unit.h IERS.cpp IERS.h Accel.cpp Accel.h EKF_GEOS3.cpp gmst.cpp gmst.h)
