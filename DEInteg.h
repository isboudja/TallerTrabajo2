//
// Created by isboudja on 15/05/2024.
//

#ifndef TALLERTRABAJO2_DEINTEG_H
#define TALLERTRABAJO2_DEINTEG_H

#include "Matrix.h"

void DEInteg(Matrix (*func)(double, Matrix &, Matrix &,Matrix &,Matrix &,Matrix &), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y,Matrix &eopdata,Matrix &Snm,Matrix &Cnm,Matrix &PC);


#endif //TALLERTRABAJO2_DEINTEG_H
