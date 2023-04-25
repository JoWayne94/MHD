/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : ExactRPSolver.cpp                          */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 23/01/2022                                                */
/*                                                                    */
/*   Programmer :                                         */
/*                                                                    */
/*   Description : Riemann Problem exact solution routines specific   */
/*                 for ideal gases.                                   */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "ExactRPSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

/**
 * @brief Constructor of the ExactRPSolver class
 */
ExactRPSolver::ExactRPSolver() {
}

/**
 * @brief Destructor of the ExactRPSolver class
 */
ExactRPSolver::~ExactRPSolver() {
}

/**
 * @brief Set variables
 *
 * @param   gamma   (Constant) adiabatic index (ratio of specific heats)
 * @param   nvar    Number of conservative variables
 */
void ExactRPSolver::SetVariables(const ADouble& wl, const ADouble& wr, double gamma, unsigned int nvar) {

    Gamma = gamma;
    nVar = nvar;

    G1 = (Gamma-1)/(2*Gamma);
    G2 = (Gamma+1)/(2*Gamma);
    G3 = 2*Gamma/(Gamma-1);
    G4 = 2/(Gamma-1);
    G5 = 2/(Gamma+1);
    G6 = (Gamma-1)/(Gamma+1);
    G7 = (Gamma-1)/2;

    rhol = wl[0];
    rhor = wr[0];
    ul = wl[1];
    ur = wr[1];
    pl = wl[4];
    pr = wr[4];

    equationofStates.SetVariables(Gamma);
    convertVars.SetVariables(nVar, Gamma, false);

    al = equationofStates.computeSoundSpeedFromEoS(rhol, pl);
    ar = equationofStates.computeSoundSpeedFromEoS(rhor, pr);
}

void ExactRPSolver::ComputeExactSolution(const double& time, const VDouble& x, const double& xd, VectArr& uExact) {

    // Check that vacuum is not generated
    double DeltaU = ur - ul;

    bool Vacuum = G4*(al+ar) <= DeltaU;
    if (Vacuum) {
        std::cerr << "The initial data is such that vacuum is generated. Program stopped.";
        std::_Exit(EXIT_FAILURE);
    }

    // Newton-Raphson method
    int iter = 0;
    int iterMax = 100;

    double Err = 1.0;
    double Tol = 1.0E-6;

    double Pold = GuessPressure();
    double Piter;

    while (Err > Tol && iter < iterMax) {

        Piter = Pold - (f_NR(Pold)/df_NR(Pold));

        Err = 2.0 * std::abs((Piter-Pold)/(Piter+Pold));

        Pold = std::max(0.0, Piter);

        iter += 1;
    }

    double Pstar = Pold;
    double ustar = 0.5*(ul+ur)+0.5*(fr_NR(Pstar)-fl_NR(Pstar));

    // Wave pattern
    int Pattern;

    if (Pstar < pl && Pstar < pr) {
        // The pattern is rarefaction-contact-rarefaction
        Pattern = 1;
    } else if (Pstar < pl && Pstar >= pr) {
        // The pattern is rarefaction-contact-shock
        Pattern = 2;
    } else if (Pstar >= pl && Pstar < pr) {
        // The pattern is shock-contact-rarefaction
        Pattern = 3;
    } else if (Pstar >= pl && Pstar >= pr) {
        // The pattern is shock-contact-shock
        Pattern = 4;
    }

    std::cout<< "Pattern = " << Pattern << " pstar = " << Pstar << " ustar = " << ustar << std::endl;

    // Sample the solution
    for (int i = 0; i < x.size(); i++) {

        double csi = (x[i] - xd)/time;

        ADouble wExact;

        wExact = Sample(Pstar, ustar, Pattern, csi);

        uExact(i, 0) = wExact;
    }
}

double ExactRPSolver::GuessPressure() const {

    double PGuess;

    // Guess pressure based on PVRS RP Solver
    int Quser = 2;
    double cup, ppv, pmin, pmax, qmax;

    cup = 0.25*(rhol+rhor)*(al+ar);
    ppv = std::max(0.0,0.5*(pl+pr)+0.5*(ul-ur)*cup);

    pmin = std::min(pl,pr);
    pmax = std::max(pl,pr);

    qmax = pmax/pmin;

    bool Cond1 = (qmax <= Quser) && ((pmin <= ppv) && (pmax >= ppv));

    if (Cond1) {
        // Select PVRS RP Solver
        PGuess = ppv;
    } else {
        if (ppv < pmin) {

            // Select two-rarefaction RP Solver
            double pq,um,ptl,ptr;

            pq = pow(pl/pr,G1);
            um = (pq*ul/al+ur/ar+G4*(pq-1.0))/(pq/al+1.0/ar);

            ptl = 1.0+G7*(ul-um)/al;
            ptr = 1.0+G7*(um-ur)/ar;

            PGuess = 0.5*(pl*pow(ptl,G3)+pr*pow(ptr,G3));

        } else {

            // Select two-shock RP solver with PVRS as estimate
            double gel,ger;

            gel = sqrt((G5/rhol)/(G6*pl+ppv));
            ger = sqrt((G5/rhor)/(G6*pr+ppv));

            PGuess = (gel*pl+ger*pr-(ur-ul))/(gel+ger);
        }
    }

    return PGuess;
}

double ExactRPSolver::fl_NR(const double& p) const {

    double fl;

    double Al = G5/rhol;
    double Bl = G6*pl;

    // Newton-Raphson function
    fl = (p > pl) ? (p-pl)*sqrt(Al/(p+Bl)) : G4*al*(pow(p/pl,G1)-1.0);

    return fl;
}

double ExactRPSolver::fr_NR(const double& p) const {

    double fr;

    double Ar = G5/rhor;
    double Br = G6*pr;

    // Newton-Raphson function
    fr = (p > pr) ? (p-pr)*sqrt(Ar/(p+Br)) : G4*ar*(pow(p/pr,G1)-1.0);

    return fr;
}

double ExactRPSolver::f_NR(const double& p) {

    double f;

    double fl = fl_NR(p);
    double fr = fr_NR(p);

    // Newton-Raphson function
    f = fl + fr + (ur - ul);

    return f;
}

double ExactRPSolver::df_NR(const double& p) const {

    double df;

    double Al = G5/rhol;
    double Ar = G5/rhor;

    double Bl = G6*pl;
    double Br = G6*pr;

    // Newton-Raphson function
    double dfL = (p > pl) ? sqrt(Al/(Bl+p))*(1-0.5*(p-pl)/(Bl+p))
                          : 1.0/(rhol*al)*pow(p/pl,-G2);
    double dfR = (p > pr) ? sqrt(Ar/(Br+p))*(1-0.5*(p-pr)/(Br+p))
                          : 1.0/(rhor*ar)*pow(p/pr,-G2);

    df = dfL + dfR;

    return df;
}

ADouble ExactRPSolver::Sample(const double& pstar, const double& ustar, const int& Pattern, const double& csi) const {

    ADouble Wexact;

    // Sampling the solution
    double rholstar, alstar, Shl, Stl, Sl;
    double rhorstar, arstar, Shr, Str, Sr;

    if (Pattern == 1 || Pattern == 2) { // r-c-r, r-c-s
        rholstar = rhol*pow(pstar/pl,1.0/Gamma); // left rarefaction wave
        alstar   = al*pow(pstar/pl,G1);
        Shl      = ul-al;
        Stl      = ustar-alstar;
    } else {
        rholstar = rhol*(pstar/pl+G6)/(1.0+G6*pstar/pl); // left shock wave
        Sl       = ul-al*sqrt(G1+G2*pstar/pl);
    }

    if (Pattern == 1 || Pattern == 3) { // r-c-r, s-c-r
        rhorstar = rhor*pow(pstar/pr,1.0/Gamma); // right rarefaction wave
        arstar   = ar*pow(pstar/pr,G1);
        Shr      = ur+ar;
        Str      = ustar+arstar;
    } else {
        rhorstar = rhor*(pstar/pr+G6)/(1.0+G6*pstar/pr); // right shock wave
        Sr       = ur+ar*sqrt(G1+G2*pstar/pr);
    }

    switch (Pattern) {
        case 1:
        {
            if (csi < Shl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[4] = pl;
            }
            else if (csi >= Shl && csi < Stl)
            {
                Wexact[0] = rhol*pow(G5+G6/al*(ul-csi),G4);
                Wexact[1] = G5*(al+G7*ul+csi);
                Wexact[4] = pl*pow(G5+G6/al*(ul-csi),G3);
            }
            else if (csi >= Stl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= ustar && csi < Str)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= Str && csi < Shr)
            {
                Wexact[0] = rhor*pow(G5-G6/ar*(ur-csi),G4);
                Wexact[1] = G5*(-ar+G7*ur+csi);
                Wexact[4] = pr*pow(G5-G6/ar*(ur-csi),G3);
            }
            else if (csi >= Shr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[4] = pr;
            }

            break;
        }

        case 2:
        {
            if (csi < Shl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[4] = pl;
            }
            else if (csi >= Shl && csi < Stl)
            {
                Wexact[0] = rhol*pow(G5+G6/al*(ul-csi),G4);
                Wexact[1] = G5*(al+G7*ul+csi);
                Wexact[4] = pl*pow(G5+G6/al*(ul-csi),G3);
            }
            else if (csi >= Stl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= ustar && csi < Sr)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= Sr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[4] = pr;
            }

            break;
        }

        case 3:
        {
            if (csi < Sl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[4] = pl;
            }
            else if (csi >= Sl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= ustar && csi < Str)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= Str && csi < Shr)
            {
                Wexact[0] = rhor*pow(G5-G6/ar*(ur-csi),G4);
                Wexact[1] = G5*(-ar+G7*ur+csi);
                Wexact[4] = pr*pow(G5-G6/ar*(ur-csi),G3);
            }
            else if (csi >= Shr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[4] = pr;
            }

            break;
        }

        case 4:
        {
            if (csi < Sl)
            {
                Wexact[0] = rhol;
                Wexact[1] = ul;
                Wexact[4] = pl;
            }
            else if (csi >= Sl && csi < ustar)
            {
                Wexact[0] = rholstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= ustar && csi < Sr)
            {
                Wexact[0] = rhorstar;
                Wexact[1] = ustar;
                Wexact[4] = pstar;
            }
            else if (csi >= Sr)
            {
                Wexact[0] = rhor;
                Wexact[1] = ur;
                Wexact[4] = pr;
            }

            default:
            {
                std::cerr << "Unrecognised wave pattern.";
                std::exit(EXIT_FAILURE);
            }
        }
    }

    return Wexact;
}