//
// Created by isboudja on 24/04/2024.
//

#ifndef UNTITLED_GLOBALS_H
#define UNTITLED_GLOBALS_H


#include "Matrix.h"


class globals {
public:
    static Matrix *eopdata;
    static Matrix *obs;
    static Matrix *Rs;
    static Matrix *Cnm;
    static Matrix *Snm;
    static Matrix *temp;
    static void eop1962();
    static void GEOS3(int nobs);
    static void GGM();
    static void DE430();

};


#endif //UNTITLED_GLOBALS_H
