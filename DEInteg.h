//
// Created by isboudja on 15/05/2024.
//

#ifndef TALLERTRABAJO2_DEINTEG_H
#define TALLERTRABAJO2_DEINTEG_H

#include "Matrix.h"

void DEInteg(Matrix (*func)(double, Matrix &), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);


#endif //TALLERTRABAJO2_DEINTEG_H
