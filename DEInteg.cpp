//
// Created by isboudja on 15/05/2024.
//

#include <limits>
#include <algorithm>
#include <iostream>
#include "DEInteg.h"
#include "Matrix.h"
#include "math.h"
#include "sign_.h"

struct DE_STATE {
    int DE_INIT = 1;
    int DE_DONE = 2;
    int DE_BADACC = 3;
    int DE_NUMSTEPS = 4;
    int DE_STIFF = 5;
    int DE_INVPARAM = 6;
};

void DEInteg(Matrix (*func)(double, Matrix &), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y) {
    double twou = 2 * std::numeric_limits<double>::epsilon();
    double fouru = 4 * std::numeric_limits<double>::epsilon();
    bool nornd;
    DE_STATE DE_STATE;

    DE_STATE = {1, 2, 3, 4, 5, 6};

    int State_ = DE_STATE.DE_INIT;
    int PermitTOUT = true;
    double told = 0;
    double ip1;
    double twoval [] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
    Matrix two(1,14,twoval,14);
    double gstr[] = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592,
                     0.00524, 0.00468};

    Matrix yy(n_eqn, 1.0);
    Matrix wt(n_eqn, 1.0);
    Matrix p(n_eqn, 1.0);
    Matrix yp(n_eqn, 1.0);
    Matrix phi(n_eqn, 17.0);
    Matrix g(14, 1.0);
    Matrix sig(14, 1.0);
    Matrix rho(14, 1.0);
    Matrix w(13, 1.0);
    Matrix alpha(13, 1.0);
    Matrix beta(13, 1.0);
    Matrix v(13, 1.0);
    Matrix psi_(13, 1.0);
    Matrix yout(n_eqn, 1.0);
    Matrix ypout(n_eqn, 1.0);
    double knew;
    double r;
    int j;
    int ifail;
    int l;
    int k;
    int i;
    int reali;
    double nsp2;
    int iq;
    double nsm2;
    double limit1;
    double limit2;
    double temp6;
    double tau;
    double  xold;
    double err;
    bool success;
    double rho2;
    double erkp1;
    Matrix F(1,1);
    if (t == tout) {
        return;
    }
    double temp1;
    double epsilon = std::max(relerr, abserr);

    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 || State_ > DE_STATE.DE_INVPARAM ||
        (State_ != DE_STATE.DE_INIT && t != told)) {
        State_ = DE_STATE.DE_INVPARAM;
        return;
    }

    double del = tout - t;
    double absdel = abs(del);
    double erkm2;
    int n;
    double nsp1;
    double erkm1;
    double hold;
    bool phase1;
    double temp2;
    double temp3;
    double temp4;
    double temp5;
    double erk;
    double absh;
    double hnew;
    double tend = t + 100.0 * del;
    bool OldPermit = true;
    double delsgn = 0.0;
    if (!PermitTOUT) {
        tend = tout;
    }
    double nostep = 0;
    double kle4 = 0;
    bool stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;
    int kold;
    double x;
    int ssj = 0;
    double h;
    bool start =false;
    int ki;
    double sum;
    double psijm1;
    double term;
    double hi;
    bool crash;
    double p5eps;
    double round;
    double gamma;
    double eta;
    int  kp1;
    int kp2;
    Matrix Dy(1,6);
    double realns;
    double im1;
    int km1;
    int km2;
    int ns;
    if (State_ == DE_STATE.DE_INIT|| !OldPermit || (delsgn * del <= 0.0)) {
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(std::max(fouru * std::abs(x), std::abs(tout - x)), tout - x);

    }

    while (true) { // Inicio del bucle de pasos

        // Si ya se ha pasado del punto de salida, interpola la solución y retorna
        if (std::abs(x - t) >= absdel) {
            for(n=1;i<=n_eqn;n++){
                yout(n, 1.0) = 0.0;
                ypout(n, 1.0) =0.0;
            }
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            hi = tout - x;
            ki = kold + 1;

            // Inicializa w[*] para calcular g[*]
            for (i = 1; i <= ki; ++i) {
                 temp1 = i;
                w(i + 1, 1) = 1.0 / temp1;
            }

            // Calcula g[*]
             term = 0.0;
            for ( j = 2; j <= ki; ++j) {
                psijm1 = psi_(j, 1);
                gamma = (hi + term) / psijm1;
               eta = hi / psijm1;
                for (i = 1; i <= ki + 1 - j; ++i) {
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                }
                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }

            // Interpola la solución yout y la derivada de la solución ypout
            for (j = 1; j <= ki; ++j) {
                    i = ki + 1 - j;
                for (int l = 0; l < phi.fil; ++l) {
                    yout = yout + g(i + 1, 1) * phi(l+1, i + 1);
                    ypout =ypout + rho(i + 1, 1) * phi(l+1, i + 1);
                }

                }



            yout =  yout*hi+y;
            y = yout;

            State_ = DE_STATE.DE_DONE; // Establece el código de

            // orno
            t = tout;                   // Establece la variable independiente
            told = t;                   // Almacena la variable independiente
            OldPermit = PermitTOUT;
            return;                     // Salida normal
        }


    if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))){
        h = tout - x;
        Dy = func(x,yy);
        for(i=1;i<=Dy.col;i++){
            yp(i,1) = Dy(1,i);
        }
        y = yy + yp * h;
        State_ = DE_STATE.DE_DONE;
        t = tout;
        told = t;
        OldPermit = PermitTOUT;
        return;
}

    h = sign_(std::min(std::abs(h), std::abs(tend - x)), h);

    for ( l = 1; l <= n_eqn; ++l) {
        wt(l,1) = releps * std::abs(yy(l,1)) + abseps;
    }

    if (abs(h) < fouru*abs(x)){
        h = sign_(fouru*abs(x),h);
        crash = true;
        return;
    }


    p5eps  = 0.5*epsilon;
    crash  = false;
    g(2,1)   = 1.0;
    g(3,1)   = 0.5;
    sig(2,1) = 1.0;

     ifail = 0;




    round = 0.0;
    for ( l=1;l<=n_eqn;l++){
            round = round + (y(l,1)*y(l,1))/(wt(l,1)*wt(l,1));
    }
            round = twou*sqrt(round);
    if (p5eps<round){
        epsilon = 2.0*round*(1.0+fouru);
        crash = true;
        return;
    }



    if(start){
        Dy = func(x,y);
        for(i=1;i<=Dy.col;i++){
            yp(i,1) = Dy(1,i);
        }
        sum = 0.0;
        for ( l=1;l<=n_eqn;l++){
        phi(l,2) = yp(l,1);
        phi(l,3) = 0.0;
        double ssj = yp(l,1)*yp(l,1);
        double ssj2 = wt(l,1)*wt(l,1);
        sum = sum + (yp(l,1)*yp(l,1))/(wt(l,1)*wt(l,1));
        }
            sum  = sqrt(sum);
        absh = abs(h);
    if (epsilon<16.0*sum*h*h){
        absh=0.25*sqrt(epsilon/sum);
    }
            h    = sign_(std::max(absh, fouru*abs(x)), h);
        hold = 0.0;
     hnew = 0.0;
     k    = 1;
     kold = 0;
    start  = false;
     phase1 = true;
         nornd  = true;
    if (p5eps<=100.0*round){
        nornd = false;
    for(l=1;l<=n_eqn;l++){
    phi(l,16)=0.0;
    }
    }
    }



    while(true){
            ssj++;
            std::cout << ssj << std::endl;
          kp1 = k+1;
         kp2 = k+2;
         km1 = k-1;
         km2 = k-2;




        if (h !=hold){
            ns=0;
        }
        if (ns<=kold){
            ns=ns+1;
        }
                nsp1 = ns+1;

        if (k>=ns){
        beta(ns+1,1) = 1.0;
         realns = ns;
        alpha(ns+1,1) = 1.0/realns;
             temp1 = h*realns;
        sig(nsp1+1,1) = 1.0;

        if (k>=nsp1){
            for (i=nsp1;i<=k;i++){

                im1   = i-1;
         temp2 = psi_(im1+1,1);
        psi_(im1+1,1) = temp1;
        beta(i+1,1)  = beta(im1+1,1)*psi_(im1+1,1)/temp2;
        temp1    = temp2 + h;
        alpha(i+1,1) = h/temp1;
        reali = i;
        sig(i+2,1) = reali*alpha(i+1,1)*sig(i+1,1);
        }
        }
        psi_(k+1,1) = temp1;



        if (ns>1){

        if (k>kold){
            temp4 = k*kp1;
        v(k+1,1) = 1.0/temp4;
         nsm2 = ns-2;
        for(j=1;j<=nsm2;j++){
                i = k-j;
        v(i+1,1) = v(i+1,1) - alpha(j+2,1)*v(i+2,1);
        }


         limit1 = kp1 - ns;
         temp5  = alpha(ns+1,1);
        for (iq=1;iq<=limit1;iq++){
        v(iq+1,1) = v(iq+1,1) - temp5*v(iq+2,1);
        w(iq+1,1) = v(iq+1,1);
        }
        g(nsp1+1,1) = w(2,1);
        }
        else{
        for (iq=1;iq<=k;iq++){

            temp3 = iq*(iq+1);
        v(iq+1,1) = 1.0/temp3;
        w(iq+1,1) = v(iq+1,1);
        }

            }
        }


         nsp2 = ns + 2;
        if (kp1>=nsp2){
            for (i=nsp2;i<=kp1;i++){
                 limit2 = kp2 - i;

       temp6  = alpha(i,1);
        for (iq=1;iq<=limit2;iq++){
        w(iq+1,1) = w(iq+1,1) - temp6*w(iq+2,1);
        }
        g(i+1,1) = w(2,1);
    }
}
        }


        //% End block 1

        // Block 2

        if (k>=nsp1){
            for (i=nsp1;i<=k;i++){
                temp1 = beta(i+1,1);
        for (l=1;i<=n_eqn;i++){
        phi(l,i+1) = temp1 * phi(l,i+1);
    }
}
}


        for (l=1;l<=n_eqn;l++){
            phi(l,kp2+1) = phi(l,kp1+1);
            phi(l,kp1+1) = 0.0;
            p(l,1)       = 0.0;
        }

        for (j=1;j<=k;j++){
                i     = kp1 - j;
         ip1   = i+1;
         temp2 = g(i+1,1);
        for (l=1;l<=n_eqn;l++){
            p(l,1)  = p(l,1) + temp2*phi(l,i+1);
            phi(l,i+1) = phi(l,i+1) + phi(l,ip1+1);
        }
        }

        if (nornd){
            p = y + p*h;
        }
        else{
            for (l=1;l<=n_eqn;l++){

                 tau = h*p(l,1) - phi(l,16);
                 p(l,1) = y(l,1) + tau;
                 phi(l,17) = (p(l,1) - y(l,1)) - tau;
        }
    }

        xold = x;
        x = x + h;
         absh = abs(h);
        Dy = func(x,p);

        for(i=1;i<=Dy.col;i++){
            yp(i,1) = Dy(1,i);
        }

         erkm2 = 0.0;
         erkm1 = 0.0;
         erk = 0.0;


        for (l=1;l<=n_eqn;l++){

                 temp3 = 1.0/wt(l,1);
                 temp4 = yp(l,1) - phi(l,1+1);
                 std::cout << yp(l,1) << std::endl;
            std::cout << phi(l,1+1) << std::endl;
        if (km2> 0){
            erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)*((phi(l,km1+1)+temp4)*temp3);
        }
        if (km2>=0){
            erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)*((phi(l,k+1)+temp4)*temp3);
    }

                erk = erk + (temp4*temp3)*(temp4*temp3);


}


        if (km2> 0){
            erkm2 = absh*sig(km1+1,1)*gstr[km2]*sqrt(erkm2);
        }
        if (km2>=0){
            erkm1 = absh*sig(k+1,1)*gstr[km1]*sqrt(erkm1);
    }



         temp5 = absh*sqrt(erk);
         err = temp5*(g(k+1,1)-g(kp1+1,1));
        erk = temp5*sig(kp1+1,1)*gstr[k];

         knew = k;

        if (km2 >0){
            if (std::max(erkm1,erkm2)<=erk){
                knew=km1;
    }
}
        if (km2==0){
            if (erkm1<=0.5*erk){
                knew=km1;
}
}
        //end block 2
        //block 3

         success = (err<=epsilon);

        if (!success){
            phase1 = false;
            x = xold;
            for (i=1;i<=k;i++){
                    temp1 = 1.0/beta(i+1,1);
                     ip1 = i+1;
            for (l=1;l<=n_eqn;l++){
            phi(l,i+1)=temp1*(phi(l,i+1)-phi(l,ip1+1));
        }
        }

            if (k>=2){
                for (i=2;i<=k;i++){
                    psi_(i,1) = psi_(i+1,1) - h;
        }
        }


            ifail = ifail+1;
            temp2 = 0.5;
            if (ifail>3){
                if (p5eps < 0.25*erk){
                    temp2 = sqrt(p5eps/erk);
                }
            }
            if (ifail>=3){
                knew = 1;
            }
                    h = temp2*h;
                    k = knew;
            if (abs(h)<fouru*abs(x)){
                crash = true;
                h = sign_(fouru*abs(x), h);
                epsilon = epsilon*2.0;
                return;
            }
        }

            if (success){
                break;
            }
        }





    kold = k;
     hold = h;



    temp1 = h*g(kp1+1,1);
    if (nornd){
        for (l=1;l<=n_eqn;l++){
            y(l,1) = p(l,1) + temp1*(yp(l,1) - phi(l,2));
        }
    }
    else{
    for (l=1;l<=n_eqn;l++){

             rho2 = temp1*(yp(l,1) - phi(l,2)) - phi(l,17);
             y(l,1) =  rho2 + p(l,1) ;
             phi(l,16) = (y(l,1) - p(l,1)) - rho2;
    }
    }
        Dy = func(x,yy);

        for(i=1;i<=Dy.col;i++){
            yp(i,1) = Dy(1,i);
        }

    for (l=1;l<=n_eqn;l++){
    phi(l,kp1+1) = yp(l,1) - phi(l,2);
    phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
    }
    for (i=1;i<=k;i++){
    for (l=1;l<=n_eqn;l++){
    phi(l,i+1) = phi(l,i+1) + phi(l,kp1+1);
}
}

     erkp1 = 0.0;


    if ( (knew==km1) || (k==12) ){
       phase1 = false;
    }

    if (phase1){
        k = kp1;
        erk = erkp1;
    }else{
    if (knew==km1){
        k = km1;
        erk = erkm1;
    }else{
    if (kp1<=ns){
        for (l=1;l<=n_eqn;l++){
            erkp1 = erkp1 + (phi(l,kp2+1)/wt(l,1))*(phi(l,kp2+1)/wt(l,1));

        }
            erkp1 = absh*gstr[kp1]*sqrt(erkp1);
    if (k>1){
        if( erkm1<=std::min(erk,erkp1)){

        k=km1; erk=erkm1;

    }else{
    if ( (erkp1<erk) && (k!=12) ){

    k=kp1;
    erk=erkp1;
    }
            }
    }else{ if (erkp1<0.5*erk){
                    k = kp1;
                        erk = erkp1;
                    }
                }
            }

        }
    }

        double ssj = two(1,k+2);
    if ( phase1 || (p5eps>=erk*two(1,k+2)) ){
        hnew = 2.0*h;

    } else{
    if (p5eps<erk){
        temp2 = k+1;
        r = p5eps/pow(erk,1.0/temp2);
        hnew = absh*std::max(0.5, std::min(0.9,r));
        hnew = sign_(std::max(hnew, fouru*abs(x)), h);
    }else{
        hnew = h;
}
}
    h = hnew;
    //std::cout << h << std::endl;
    //std::cout << two(1,k+2) << std::endl;
    //std::cout << erk << std::endl;

    if (crash){
        State_    = DE_STATE.DE_BADACC;
        relerr    = epsilon*releps;
        abserr    = epsilon*abseps;
        y         = yy;
        t         = x;
        told      = t;
        OldPermit = true;
        return;
    }

    nostep = nostep+1;


    kle4 = kle4+1;
    if (kold>  4){
        kle4 = 0;
    }
    if (kle4>=50){
        stiff = true;
    }
}



}









