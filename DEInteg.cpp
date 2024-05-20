//
// Created by isboudja on 15/05/2024.
//

#include <limits>
#include <algorithm>
#include "DEInteg.h"
#include "Matrix.h"
#include "math.h"
#include "sign_.h"
/*
struct DE_STATE {
    int DE_INIT = 1;
    int DE_DONE = 2;
    int DE_BADACC = 3;
    int DE_NUMSTEPS = 4;
    int DE_STIFF = 5;
    int DE_INVPARAM = 6;
};

Matrix DEInteg(double (*func)(double, Matrix), double t, double tout, double relerr, double abserr, int n_eqn, Matrix y) {
    double twou = 2 * std::numeric_limits<double>::epsilon();
    double fouru = 4 * std::numeric_limits<double>::epsilon();

    DE_STATE DE_STATE;

    DE_STATE = {1, 2, 3, 4, 5, 6};

    int State_ = DE_STATE.DE_INIT;
    int PermitTOUT = true;
    double told = 0;

    double two[] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
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

    if (t == tout) {

    }

    double epsilon = std::max(relerr, abserr);

    if (relerr < 0.0 || abserr < 0.0 || epsilon <= 0.0 || State_ > DE_STATE.DE_INVPARAM ||
        (State_ != DE_STATE.DE_INIT && t != told)) {
        State_ = DE_STATE.DE_INVPARAM;

    }

    double del = tout - t;
    double absdel = abs(del);

    double tend = t + 100.0 * del;
    bool OldPermit = false;
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
    double h;
    if (State_ == DE_STATE.DE_INIT || !OldPermit || (delsgn * del <= 0.0)) {
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        bool start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(std::max(fouru * std::abs(x), std::abs(tout - x)), tout - x);
    }

    while (true) { // Inicio del bucle de pasos

        // Si ya se ha pasado del punto de salida, interpola la solución y retorna
        if (std::abs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1.0);
            Matrix ypout(n_eqn, 1.0);
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            double hi = tout - x;
            int ki = kold + 1;

            // Inicializa w[*] para calcular g[*]
            for (int i = 1; i <= ki; ++i) {
                double temp1 = i;
                w(i + 1, 1) = 1.0 / temp1;
            }

            // Calcula g[*]
            double term = 0.0;
            for (int j = 2; j <= ki; ++j) {
                double psijm1 = psi_(j, 1);
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 1; i <= ki + 1 - j; ++i) {
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                }
                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }

            // Interpola la solución yout y la derivada de la solución ypout
            for (int j = 1; j <= ki; ++j) {
                int i = ki + 1 - j;
                for (int l = 0; l < n_eqn; ++l) {
                    yout(l, 1) += g(i + 1, 1) * phi(l, i + 1);
                    ypout(l, 1) += rho(i + 1, 1) * phi(l, i + 1);
                }
            }

            for (int l = 0; l < n_eqn; ++l) {
                yout(l, 1) = y(l, 1) + hi * yout(l, 1);
                y(l, 1) = yout(l, 1);
            }
            State_ = DE_STATE.DE_DONE; // Establece el código de retorno
            t = tout;                   // Establece la variable independiente
            told = t;                   // Almacena la variable independiente
            OldPermit = PermitTOUT;
            break;                     // Salida normal
        }
    }

    if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))){
        h = tout - x;
        yp = func(x, yy);
        y = yy + yp * h;
        State_ = DE_STATE.DE_DONE;
        t = tout;
        told = t;
        OldPermit = PermitTOUT;
        return;
}



    /*Y =

            5542555.8942746
    3213514.83814162
    3990892.92789062
    5394.06894044333
    -2365.21290573999
    -7061.84481373379




    return y;
}*/