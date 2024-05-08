//
// Created by itsfr on 02/05/2024.
//

#include <cmath>
#include <iostream>
#include "IERS.h"
#include "SAT_Const.h"
Matrix IERS(Matrix& eop, double Mjd_UTC, char interp){

double mfme, fixf;
Matrix Result(1,9);
Matrix preeop(1,13);
Matrix nexteop(1,13);
int i;
if (interp == 'l') {
int mjd = floor(Mjd_UTC);
    for (i = 0; i < 13; i++) {
        if(eop(4,i+1)==mjd){
            mjd = i+1;
            break;

        }
    }
for (i = 0; i < 13; i++) {
preeop(1,i+1) = eop(i+1,mjd);
nexteop(1,i+1) = eop(i+1,mjd+1);
}
mfme = 1440. * (Mjd_UTC-floor(Mjd_UTC));
fixf = mfme / 1440.;

    Result(1,1) = preeop(1,5) + (nexteop(1,5) - preeop(1,5)) * fixf;
    Result(1,2) = preeop(1,6) + (nexteop(1,6) - preeop(1,6)) * fixf;
    Result(1,3) = preeop(1,7) + (nexteop(1,7) - preeop(1,7)) * fixf;
    Result(1,4)= preeop(1,8) + (nexteop(1,8) - preeop(1,8))* fixf;
    Result(1,5) = preeop(1,9) + (nexteop(1,9) - preeop(1,9)) * fixf;
    Result(1,6) = preeop(1,10) + (nexteop(1,10) - preeop(1,10)) * fixf;
    Result(1,7)= preeop(1,11) + (nexteop(1,11) - preeop(1,11)) * fixf;
    Result(1,8) = preeop(1,12) + (nexteop(1,12) - preeop(1,12)) * fixf;
    Result(1,9) =  preeop(1,13);

    Result(1,1) /= Constants::Arcs*1.;
    Result(1,2) /= Constants::Arcs*1.;
    Result(1,5) /= Constants::Arcs*1.;
    Result(1,6) /= Constants::Arcs*1.;
    Result(1,7) /= Constants::Arcs*1.;
    Result(1,8) /= Constants::Arcs*1.;

    preeop.print();
} else if (interp == 'n') {
int mjd = floor(Mjd_UTC);
    for (i = 0; i < 13; i++) {
        if(eop(4,i+1)==mjd){
            mjd = i+1;
            break;

        }
    }
for (i = 0; i < 13; i++) {
    preeop(1,i+1) = eop(mjd,i+1);
}

    Result(1,1)  = preeop(1,5) / Constants::Arcs;
    Result(1,2) = preeop(1,6) / Constants::Arcs;
    Result(1,3)= preeop(1,7);
    Result(1,4) = preeop(1,8);
    Result(1,5)= preeop(1,9) / Constants::Arcs;
    Result(1,6) = preeop(1,10) / Constants::Arcs;
    Result(1,7) = preeop(1,11) / Constants::Arcs;
    Result(1,8) = preeop(1,12) / Constants::Arcs;
    Result(1,9) = preeop(1,13);

}

    return Result;

}