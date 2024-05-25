//
// Created by isboudja on 22/05/2024.
//

#include <cstdio>
#include "AccelHarmonic.h"
#include "SAT_Const.h"
#include "Legendre.h"

Matrix AccelHarmonic(Matrix &r,Matrix &E,double n_max,double m_max) {

    FILE *fid = fopen("../texts/GGM03S.txt","r");
    Matrix Cnm(181,181);
    Matrix Snm(181,181);
    Matrix temp(1,6);
    Matrix a_bf(1,3);
    double c;
    if(fid==nullptr){
        printf("error");
        exit(EXIT_FAILURE);
    }
    int n;
    int m;
    for (n=0;n<=180;n++){
        for (m=0;m<=n;m++) {
            fscanf(fid, "%d%d%lf%lf%lf%lf",&c,&c,&temp(1,3),&temp(1,4),&temp(1,5),&temp(1,6));
            Cnm(n+1,m+1) = temp(1,3);
            Snm(n+1,m+1) = temp(1,4);
        }
    }
    Constants constants;

    double r_ref = 6378.1363e3;
    double gm = 398600.4415e9;
    Matrix r_bf =  E*r;
    Matrix r2(1,3);

    r_bf.print();
    if(r_bf.col<3){
        for( n=1;n<=3;n++){
            r2(1,n) = r_bf(n,1);
        }
    }
    r2.print();
    double d = r2.norm();
    double latgc = asin(r_bf(3,1) / d);
    double lon = atan2(r_bf(2,1), r_bf(1,1));

     Matrix Leg = Legendre(n_max, m_max, latgc);
     Leg.print();
     Matrix pnm(1,(n_max+1)*(m_max+1));
     Matrix dpnm(1,(n_max+1)*(m_max+1));
     for(int i=1;i<=(n_max+1)*(m_max+1);i++){
         pnm(1,i) = Leg(1,i);
     }
     pnm.print();
     int j = 1;
    for(int i=(n_max+1)*(m_max+1)+1;i<=((n_max+1)*(m_max+1))*2;i++){
        dpnm(1,j) = Leg(1,i);
        j++;
    }
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;
    double dUdr = 0.0, dUdlatgc = 0.0, dUdlon = 0.0;
    double b1,b2,b3;
    j = 1;
    for(double n = 0.; n <= n_max; n++) {
         b1 = (-gm / pow(d, 2)) * pow((r_ref / d), n) * (n + 1);
         b2 = (gm / d) * pow((r_ref / d), n);
         b3 = (gm / d) * pow((r_ref / d), n);

        for(double m = 0.; m <= m_max; m++) {
            double pnm_val = pnm(1, j); // Almacenar el valor para evitar recalculaciones
            double dpnm_val = dpnm(1, j); // Almacenar el valor para evitar recalculaciones

            double Cnm_val = Cnm(n + 1, m + 1); // Almacenar el valor para evitar recalculaciones
            double Snm_val = Snm(n + 1, m + 1); // Almacenar el valor para evitar recalculaciones

            q1 += pnm_val * (Cnm_val * cos(m * lon) + Snm_val * sin(m * lon));
            q2 += dpnm_val * (Cnm_val * cos(m * lon) + Snm_val * sin(m * lon));
            q3 += m * pnm_val * (Snm_val * cos(m * lon) - Cnm_val * sin(m * lon));
            j += 1;

        }
        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;

        // Restablecer q1, q2, y q3 a cero para la próxima iteración de n
        q1 = 0.0;
        q2 = 0.0;
        q3 = 0.0;
    }

    double r2xy = pow(r_bf(1, 1),2) + pow(r_bf(2, 1),2);
    double ax = (1. / d * dUdr - r_bf(3,1) / (pow(d,2) * sqrt(r2xy)) * dUdlatgc) * r_bf(1,1) - (1. / r2xy * dUdlon) * r_bf(2,1);
    double ay = (1. / d * dUdr - r_bf(3,1) / (pow(d,2) * sqrt(r2xy)) * dUdlatgc) * r_bf(2,1) + (1. / r2xy * dUdlon) * r_bf(1,1);
    double az = 1./ d * dUdr * r_bf(3,1) + sqrt(r2xy) / pow(d,2) * dUdlatgc;


    a_bf(1,1) = ax;
    a_bf(1,2) = ay;
    a_bf(1,3) = az;



    a_bf.print();
    E.print();
    Matrix E2 = E.transpose();
    E2.print();

    Matrix a = a_bf*E;
    a.print();
    return a;
}



