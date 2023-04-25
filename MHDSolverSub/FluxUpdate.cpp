/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : FluxUpdate.cpp                             */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 21/01/2022                                                */
/*                                                                    */
/*   Programmer :                                         */
/*                                                                    */
/*   Description : Compute stable time step size for each iteration.  */
/*                 Define class member functions for the FluxUpdate   */
/*                 class and update the solution to t ^ {n + 1}       */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "FluxUpdate.h"
#include <iostream>
#include <cmath>

using namespace std;

/**
 * @brief Constructor of the FluxUpdate class
 */
FluxUpdate::FluxUpdate() {
}

/**
 * @brief Destructor of the FluxUpdate class
 *
 * 	      Clean up memory
 */
FluxUpdate::~FluxUpdate() {
    /// Deallocate memory of member vectors
    delete Ql;
    delete Qr;
}

/**
 * @brief Set variables
 */
void FluxUpdate::SetVariables(unsigned int nx, unsigned int ny, double cfl, unsigned int nxGhost, unsigned int nyGhost,
                              unsigned int nvar, double finalt, double gamma, bool isdivclean) {
    Nx = nx;
    Ny = ny;
    Cfl = cfl;
    NxGhost = nxGhost;
    NyGhost = nyGhost;
    nVar = nvar;
    T = finalt;
    Gamma = gamma;

    numericalFlux.SetVariables(nVar, Gamma, isdivclean);
    convertVars.SetVariables(nVar, Gamma, isdivclean);
    equationofStates.SetVariables(Gamma);
}

/**
 * @brief Initialise Ql and Qr vectors
 *
 *	      Set conservative vectors for left and right states.
 *	      Can be at half time step if slope-limiting is employed
 */
void FluxUpdate::Initialise() {

    Ql = new VectArr;
    Ql->SetSize(Nx + NxGhost, Ny + NyGhost);
    Qr = new VectArr;
    Qr->SetSize(Nx + NxGhost, Ny + NyGhost);

}

/**
 * @brief Compute time step (dt) routine
 *
 *        Toro's book pg 221 section 6.3.2
 *
 * @param   ubc      Vector of conservative variables including fictitious cells
 * @param   dx       Spatial discretisation length in x-direction
 * @param   dy       Spatial discretisation length in y-direction
 * @param   time     Current time
 * @param   dt       Time step size
 */
void FluxUpdate::ComputeDt(const VectArr& ubc, const double& deltax, const double& deltay, const double& time, double& deltat, double& Smax) {

    double Slx, Srx, Sly, Sry; // wave speeds at left right bottom top
    ADouble qlx, qrx; // Left and right conservative states of the RP
    ADouble qly, qry; // Bottom and top conservative states of the RP
    ADouble w;

    for (int i = 0; i < Nx + 1; i++) {
        for (int j = 0; j < Ny + 1; j++) {

            qlx = ubc(i + NxGhost - 1, j + NyGhost - 1);
            qrx = ubc(i + NxGhost, j + NyGhost - 1);

            qly = ubc(i + NxGhost - 1, j + NyGhost - 1);
            qry = ubc(i + NxGhost - 1, j + NyGhost);

//            w = convertVars.conservativeToprimitive((ADouble) ubc(i, j));
//
//            double rho = w[0];
//            double ux = w[1];
//            double uy = w[2];
//            double uz = w[3];
//            double p = w[4];
//            double Bx = w[5];
//            double By = w[6];
//            double Bz = w[7];
//
//            double a = equationofStates.computeSoundSpeedFromEoS((double) rho, (double) p);
//
//            double ca = sqrt(Bx*Bx + By*By + Bz*Bz)/sqrt(rho); // Strictly not the Alfvén wave speed but for convenience defined as such
//            double cfx, cfy;
//
//            cfx = sqrt(0.5 * (a*a + ca*ca + sqrt((a*a + ca*ca)*(a*a + ca*ca) - 4 * (a*a * Bx*Bx/rho))));
//            cfy = sqrt(0.5 * (a*a + ca*ca + sqrt((a*a + ca*ca)*(a*a + ca*ca) - 4 * (a*a * By*By/rho))));
//
//            Smax = std::max(Smax, std::max(fabs(ux) + cfx, fabs(uy) + cfy));

             numericalFlux.WaveEstimates(qlx, qrx, Slx, Srx, 1); // get Sl and Sr
             numericalFlux.WaveEstimates(qly, qry, Sly, Sry, 2);

             Smax = std::max(Smax, std::max({fabs(Slx), fabs(Sly), fabs(Srx), fabs(Sry)})); // eq (6.19)
        }
    }

    deltat = std::min(Cfl * std::min<double>(deltax, deltay) / Smax, T - time); // eq (6.17), make sure dt addition does not exceed final time output
}

/**
 * @brief Half time step local evolution routine
 *
 *        CCM1 Lecture 6
 *
 * @param   ubc      Vector of conservative variables including fictitious cells
 * @param   dx       Spatial discretisation length in x-direction
 * @param   dy       Spatial discretisation length in y-direction
 * @param   time     Current time
 * @param   dt       Time step size
 * @param   k, l     Current iteration number
 */
void FluxUpdate::HalfTimeStepUpdate(const VectArr& ubc, const double& deltax, const double& deltat, ADouble& qhtsl, ADouble& qhtsr, int k, int l, int dim, const double& ch) {

    double omega = 0.0;
    ADouble qll, ql, qr;
    ADouble r, etaR, eta, deltaminus, deltaplus, delta;
    ADouble Qbarl, Qbarr, Fl, Fr;

    if (dim == 1) {
        qll = ubc(k + NxGhost - 2, l + NyGhost - 1); // i-1, j
        ql = ubc(k + NxGhost - 1, l + NyGhost - 1);  // i,   j
        qr = ubc(k + NxGhost, l + NyGhost - 1);      // i+1, j
    } else if (dim == 2) {
        qll = ubc(k + NxGhost - 1, l + NyGhost - 2); // i, j-1
        ql = ubc(k + NxGhost - 1, l + NyGhost - 1);  // i, j
        qr = ubc(k + NxGhost - 1, l + NyGhost);      // i, j+1
    }

    /// Minbee
    for (int i = 0; i < nVar; i++) {
        r[i] = (ql[i] - qll[i]) / (qr[i] - ql[i]);
        etaR[i] = 2 / (1 + r[i]);
        if (r[i] <= 0) { eta[i] = 0; }
        else if (r[i] > 0 && r[i] <= 1) { eta[i] = r[i]; }
        else if (r[i] > 1) { eta[i] = std::min<double>(1, etaR[i]); }

        deltaminus[i] = ql[i] - qll[i];
        deltaplus[i] = qr[i] - ql[i];
        delta[i] = 0.5 * (1 + omega) * deltaminus[i] + 0.5 * (1 - omega) * deltaplus[i];
        Qbarl[i] = ql[i] - 0.5 * eta[i] * delta[i];
        Qbarr[i] = ql[i] + 0.5 * eta[i] * delta[i];
    }

    /// Compute the MHD flux
    Fl = convertVars.MhdFlux((ADouble) Qbarl, dim, ch);
    Fr = convertVars.MhdFlux((ADouble) Qbarr, dim, ch);

    for (int i = 0; i < nVar; i++) {
        qhtsl[i] = Qbarl[i] - 0.5 * (deltat/deltax) * (Fr[i] - Fl[i]);
        qhtsr[i] = Qbarr[i] - 0.5 * (deltat/deltax) * (Fr[i] - Fl[i]);
    }

//    if (dim == 1) {
//        qhtsl[5] = Qbarl[5];
//        qhtsr[5] = Qbarr[5];
//        qhtsl[8] = Qbarl[8];
//        qhtsr[8] = Qbarr[8];
//    } else if (dim == 2) {
//        qhtsl[6] = Qbarl[6];
//        qhtsr[6] = Qbarr[6];
//        qhtsl[8] = Qbarl[8];
//        qhtsr[8] = Qbarr[8];
//    }
}

/**
 * @brief Update routine
 *
 *        Toro's book pg 175 section 5.3.1
 *
 * @param   ubcOld      Old vector of conservative variables including fictitious cells
 * @param   unew        New vector of conservative variables
 * @param   dx          Spatial discretisation length in x-direction
 * @param   dy          Spatial discretisation length in y-direction
 * @param   dt          Time step size
 */
void FluxUpdate::UpdatewithFluxes(const VectArr& ubcOld, VectArr& unew, const double& deltax, const double& deltat, int dim, const double& ch, PairDouble& timerPair) {

    VectArr NumericalFluxes;
    NumericalFluxes.SetSize(Nx + 1, Ny + 1);

    start = clock();

    /// Compute the numerical flux
    for (int i = 0; i < Nx + NxGhost; i++) {
        for (int j = 0; j < Ny + NyGhost; j++) {

            /// If not slope-limiting
            // Ql = ubcOld[i + NxGhost - 1];
            // Qr = ubcOld[i + NxGhost];

            /// Compute Half time step evolution (slope-limiting)
            HalfTimeStepUpdate(ubcOld, deltax, deltat, (*Ql)(i, j), (*Qr)(i, j), i, j, dim, ch); // Output Ql and Qr
        }
    }

    end = clock();
    elapsed_timehts = (end - start)/(double)CLOCKS_PER_SEC ;
    timerPair.first = elapsed_timehts;

    /// Implement and replace here other numerical fluxes, for example HLL, HLLC, FORCE or TV. Outputs NumericalFluxes[i]
    // numericalFlux->FORCEFlux((*Qr)[i][j+1], (*Ql)[i+1][j+1], NumericalFluxes[i][j], deltax, deltat, dim);
    // numericalFlux->FORCEFlux((*Qr)[i+1][j], (*Ql)[i+1][j+1], NumericalFluxes[i][j], deltax, deltat, dim);
    ADouble Fl, Fr;

    start = clock();

    if (dim == 1) {

        for (int i = 0; i < Nx + 1; i++) {
            for (int j = 0; j < Ny; j++) {
                numericalFlux.HLLCFlux((*Qr)(i, j+1), (*Ql)(i+1, j+1), NumericalFluxes(i, j), dim, ch);
                // numericalFlux.FORCEFlux((*Qr)(i, j+1), (*Ql)(i+1, j+1), NumericalFluxes(i, j), deltax, deltat, dim, ch);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Fl = NumericalFluxes(i, j);
                Fr = NumericalFluxes(i+1, j);
                /** Conservative method, Eq (5.42) **/
                for (int n = 0; n < nVar; n++) {
                    unew(i, j, n) = ubcOld(i + NxGhost, j + NyGhost, n) + deltat/deltax * (Fl[n] - Fr[n]);
                }
            }
        }

    } else if (dim == 2) {

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny + 1; j++) {
                numericalFlux.HLLCFlux((*Qr)(i+1, j), (*Ql)(i+1, j+1), NumericalFluxes(i, j), dim, ch);
                // numericalFlux.FORCEFlux((*Qr)(i+1, j), (*Ql)(i+1, j+1), NumericalFluxes(i, j), deltax, deltat, dim, ch);
            }
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Fl = NumericalFluxes(i, j);
                Fr = NumericalFluxes(i, j+1);
                /** Conservative method, Eq (5.42) **/
                for (int n = 0; n < nVar; n++) {
                    unew(i, j, n) = ubcOld(i + NxGhost, j + NyGhost, n) + deltat/deltax * (Fl[n] - Fr[n]);
                }
            }
        }
    }

    end = clock();
    elapsed_timeflux = (end - start)/(double)CLOCKS_PER_SEC ;
    timerPair.second = elapsed_timeflux;

}

void FluxUpdate::OperatorSplitting(const VectArr& ubcOld, VectArr& unew, const double& deltat, const double& ch) const {

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            unew(i, j, 8) = ubcOld(i + NxGhost, j + NyGhost, 8) * exp(- deltat * (ch / 0.18)); // Dedner et al.
        }
    }

}

// VVDouble temp(Nx + 1, VDouble(Ny + 1));
// cs = equationofStates->computeSoundSpeedFromEoS((double) ubc[i + NxGhost - 1][j + NyGhost - 1][0], (double) ubc[i + NxGhost - 1][j + NyGhost - 1][3]);
// temp[i][j] = sqrt(ubc[i + NxGhost - 1][j + NyGhost - 1][1]*ubc[i + NxGhost - 1][j + NyGhost - 1][1] + ubc[i + NxGhost - 1][j + NyGhost - 1][2]*ubc[i + NxGhost - 1][j + NyGhost - 1][2] + ubc[i + NxGhost - 1][j + NyGhost - 1][2]*ubc[i + NxGhost - 1][j + NyGhost - 1][3]) + cs;
// Smax = std::max(Smax, temp[i][j]); // eq (6.20)
