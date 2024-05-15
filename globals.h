//
// Created by isboudja on 24/04/2024.
//

#ifndef UNTITLED_GLOBALS_H
#define UNTITLED_GLOBALS_H


#include "Matrix.h"


class globals {
public:
    static Matrix *matrix;
    static Matrix *matrix2;
    static Matrix *Cnm;
    static Matrix *Snm;
    static Matrix *temp;
    static void eop1962(int c);
    static void GEOS3(int c);
    static void GGM(int c);
    static void DE430();

};


#endif //UNTITLED_GLOBALS_H
