/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : NumericalFlux.cpp                          */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 22/01/2022                                                */
/*                                                                    */
/*   Programmer :                                         */
/*                                                                    */
/*   Description : Various centred and Riemann schemes to             */
/*                 compute the numerical fluxes at i+1/2 and i-1/2    */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include "NumericalFlux.h"

using namespace std;

/**
 * @brief Constructor of the NumericalFlux class
 */
NumericalFlux::NumericalFlux() {
}

/**
 * @brief Destructor of the NumericalFlux class
 */
NumericalFlux::~NumericalFlux() {
}

/**
 * @brief Set variables
 *
 * @param   nvar     Number of conservative variables
 * @param   gamma    (Constant) adiabatic index (ratio of specific heats)
 */
void NumericalFlux::SetVariables(unsigned int nvar, double gamma, bool isdivclean) {

    nVar = nvar;
    Gamma = gamma;
    isDivClean = isdivclean;

    convertVars.SetVariables(nVar, Gamma, isdivclean);
    equationofStates.SetVariables(Gamma);

}

/**
 * @brief Estimate wave speeds S_L and S_R in x and y direction
 *
 *        Toro's book pg 327 section 10.5
 *
 * @param   Ql      Left vector of conservative variables
 * @param   Qr      Right vector of conservative variables
 * @param   Sl      Left wave speed
 * @param   Sr      Right wave speed
 */
void NumericalFlux::WaveEstimates(const ADouble& Ql, const ADouble& Qr, double& Sleft, double& Sright, int dim) {

    ADouble wl, wr;
    wl = convertVars.conservativeToprimitive((ADouble) Ql);
    wr = convertVars.conservativeToprimitive((ADouble) Qr);

    double rhol = wl[0];
    double rhor = wr[0];
    double uxl = wl[1];
    double uxr = wr[1];
    double uyl = wl[2];
    double uyr = wr[2];
    double uzl = wl[3];
    double uzr = wr[3];
    double pl = wl[4];
    double pr = wr[4];
    double Bxl = wl[5];
    double Bxr = wr[5];
    double Byl = wl[6];
    double Byr = wr[6];
    double Bzl = wl[7];
    double Bzr = wr[7];

    /** Ideal Gas EoS **/
    double al = equationofStates.computeSoundSpeedFromEoS((double) rhol, (double) pl);
    double ar = equationofStates.computeSoundSpeedFromEoS((double) rhor, (double) pr);

    double cal = sqrt(Bxl*Bxl + Byl*Byl + Bzl*Bzl)/sqrt(rhol); // Strictly not the Alfvén wave speed but for convenience defined as such
    double car = sqrt(Bxr*Bxr + Byr*Byr + Bzr*Bzr)/sqrt(rhor);
    double cfl, cfr;

    if (dim == 1) {

        cfl = sqrt(0.5 * (al*al + cal*cal + sqrt((al*al + cal*cal)*(al*al + cal*cal) - 4 * (al*al * Bxl*Bxl/rhol))));
        cfr = sqrt(0.5 * (ar*ar + car*car + sqrt((ar*ar + car*car)*(ar*ar + car*car) - 4 * (ar*ar * Bxr*Bxr/rhor))));

        Sleft = min(uxl, uxr) - max(cfl, cfr);
        Sright = max(uxl, uxr) + max(cfl, cfr);

    } else if (dim == 2) {

        cfl = sqrt(0.5 * (al*al + cal*cal + sqrt((al*al + cal*cal)*(al*al + cal*cal) - 4 * (al*al * Byl*Byl/rhol))));
        cfr = sqrt(0.5 * (ar*ar + car*car + sqrt((ar*ar + car*car)*(ar*ar + car*car) - 4 * (ar*ar * Byr*Byr/rhor))));

        Sleft = min(uyl, uyr) - max(cfl, cfr);
        Sright = max(uyl, uyr) + max(cfl, cfr);

    }

    /// More restrictive wave speed definitions
//     Sleft = min(sqrt(uxl*uxl + uyl*uyl + uzl*uzl), sqrt(uxr*uxr + uyr*uyr + uzr*uzr)) - max(cfl, cfr);
//     Sright = max(sqrt(uxl*uxl + uyl*uyl + uzl*uzl), sqrt(uxr*uxr + uyr*uyr + uzr*uzr)) + max(cfl, cfr);

    // Add here different Wave speed estimates ...
}

/**
 * @brief HLL Approx Riemann Solver
 *
 *        Toro's book pg 320 section 10.3
 *
 * @param   Ql        Left vector of conservative variables
 * @param   Qr        Right vector of conservative variables
 * @param   Fhll      HLL intercell flux for the approximate Godunov method
 */
void NumericalFlux::HLLFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fhll, int dim, const double& ch) {

    WaveEstimates(Ql, Qr, Sl, Sr, dim); // calculate Sl and Sr

    ADouble Fl, Fr;

    Fl = convertVars.MhdFlux((ADouble) Ql, dim, ch);
    Fr = convertVars.MhdFlux((ADouble) Qr, dim, ch);

    /** Eq (10.21) **/
    if (Sl >= 0) {
        Fhll = Fl;
        return;
    }

    if (Sr <= 0) {
        Fhll = Fr;
        return;
    }

    for (int n = 0; n < nVar; n++) {
        Fhll[n] = (Sr * Fl[n] - Sl * Fr[n] + Sl * Sr * (Qr[n] - Ql[n]))/(Sr - Sl);
    }
}

/**
 * @brief HLLC Approx Riemann Solver
 *
 *        Toro's book pg 322 section 10.4
 *
 * @param   Ql         Left vector of conservative variables
 * @param   Qr         Right vector of conservative variables
 * @param   Fhllc      HLLC intercell flux for the approximate Godunov method
 */
void NumericalFlux::HLLCFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fhllc, int dim, const double& ch) {

    double Sstar;
    ADouble wl, wr, Ustarl, Ustarr, Wstarl, Wstarr;

    wl = convertVars.conservativeToprimitive((ADouble) Ql);
    wr = convertVars.conservativeToprimitive((ADouble) Qr);

    double rhol = wl[0];
    double rhor = wr[0];
    double uxl = wl[1];
    double uxr = wr[1];
    double uyl = wl[2];
    double uyr = wr[2];
    double uzl = wl[3];
    double uzr = wr[3];
    double pl = wl[4];
    double pr = wr[4];
    double Bxl = wl[5];
    double Bxr = wr[5];
    double Byl = wl[6];
    double Byr = wr[6];
    double Bzl = wl[7];
    double Bzr = wr[7];
    double psil, psir;

    ADouble Qlnew, Qrnew;
    double Btildel, Btilder, psitildel, psitilder;
    Qlnew = Ql;
    Qrnew = Qr;

    if (isDivClean) {
        psil = wl[8];
        psir = wr[8];

        if (dim == 1) {
            Btildel = 0.5 * (Bxl + Bxr) - (0.5/ch) * (psir - psil);
            Btilder = 0.5 * (Bxl + Bxr) - (0.5/ch) * (psir - psil);
            psitildel = 0.5 * (psil + psir) - (ch/2.0) * (Bxr - Bxl);
            psitilder = 0.5 * (psil + psir) - (ch/2.0) * (Bxr - Bxl);
            Qlnew[5] = Btildel;
            Qrnew[5] = Btilder;
        } else if (dim == 2) {
            Btildel = 0.5 * (Byl + Byr) - (0.5/ch) * (psir - psil);
            Btilder = 0.5 * (Byl + Byr) - (0.5/ch) * (psir - psil);
            psitildel = 0.5 * (psil + psir) - (ch/2.0) * (Byr - Byl);
            psitilder = 0.5 * (psil + psir) - (ch/2.0) * (Byr - Byl);
            Qlnew[6] = Btildel;
            Qrnew[6] = Btilder;
        }

        Qlnew[8] = psitildel;
        Qrnew[8] = psitilder;
    } else {
        if (dim == 1) {
            Btildel = Bxl;
            Btilder = Bxr;
        } else if (dim == 2) {
            Btildel = Byl;
            Btilder = Byr;
        }
    }

    WaveEstimates(Qlnew, Qrnew, Sl, Sr, dim); // calculate Sl and Sr

    ADouble Fl, Fr, Qstarl, Qstarr, Fstarl, Fstarr;

    Fl = convertVars.MhdFlux((ADouble) Qlnew, dim, ch);
    Fr = convertVars.MhdFlux((ADouble) Qrnew, dim, ch);

    /** Eq (10.26) **/
    if (Sl >= 0) {
        Fhllc = Fl;
        return;
    }

    if (Sr <= 0) {
        Fhllc = Fr;
        return;
    }

    for (int n = 0; n < nVar; n++) {
        Ustarl[n] = (Sr * Qrnew[n] - Sl * Qlnew[n] - (Fr[n] - Fl[n]))/(Sr - Sl);
        Ustarr[n] = Ustarl[n];
    }

    Wstarl = convertVars.conservativeToprimitive(Ustarl);
    Wstarr = convertVars.conservativeToprimitive(Ustarr);

    double BhllvhllL = Wstarl[5]*Wstarl[1] + Wstarl[6]*Wstarl[2] + Wstarl[7]*Wstarl[3];
    double BhllvhllR = Wstarr[5]*Wstarr[1] + Wstarr[6]*Wstarr[2] + Wstarr[7]*Wstarr[3];

    double pTLstar, pTRstar, rhoLstar, rhoRstar, rhouxLstar, rhouxRstar, rhouyLstar, rhouyRstar, rhouzLstar, rhouzRstar, ELstar, ERstar;

    if (dim == 1) {

        double pTl = pl + 0.5 * (Btildel*Btildel + Byl*Byl + Bzl*Bzl);
        double pTr = pr + 0.5 * (Btilder*Btilder + Byr*Byr + Bzr*Bzr);

        double BdotvL = Btildel*uxl + Byl*uyl + Bzl*uzl;
        double BdotvR = Btilder*uxr + Byr*uyr + Bzr*uzr;

        Sstar = (pTr - pTl - Btildel*Btildel + Btilder*Btilder + rhol * uxl * (Sl - uxl) - rhor * uxr * (Sr - uxr)) / (rhol * (Sl - uxl) - rhor * (Sr - uxr)); // eq (10.37)

        pTLstar = rhol*(Sl - uxl)*(Sstar - uxl) + pTl - Btildel*Btildel + Wstarl[5]*Wstarl[5];
        pTRstar = rhor*(Sr - uxr)*(Sstar - uxr) + pTr - Btilder*Btilder + Wstarr[5]*Wstarr[5];

        rhoLstar = rhol * (Sl - uxl)/(Sl - Sstar);
        rhoRstar = rhor * (Sr - uxr)/(Sr - Sstar);

        rhouxLstar = rhoLstar * Sstar;
        rhouxRstar = rhoRstar * Sstar;

        rhouyLstar = rhol*uyl * (Sl - uxl)/(Sl - Sstar) - (Wstarl[5]*Wstarl[6] - Btildel*Byl)/(Sl - Sstar);
        rhouyRstar = rhor*uyr * (Sr - uxr)/(Sr - Sstar) - (Wstarr[5]*Wstarr[6] - Btilder*Byr)/(Sr - Sstar);

        rhouzLstar = rhol*uzl * (Sl - uxl)/(Sl - Sstar) - (Wstarl[5]*Wstarl[7] - Btildel*Bzl)/(Sl - Sstar);
        rhouzRstar = rhor*uzr * (Sr - uxr)/(Sr - Sstar) - (Wstarr[5]*Wstarr[7] - Btilder*Bzr)/(Sr - Sstar);

        ELstar = Ql[4] * (Sl - uxl)/(Sl - Sstar) + (pTLstar*Sstar - pTl*uxl - (Wstarl[5] * BhllvhllL - Btildel * BdotvL))/(Sl - Sstar);
        ERstar = Qr[4] * (Sr - uxr)/(Sr - Sstar) + (pTRstar*Sstar - pTr*uxr - (Wstarr[5] * BhllvhllR - Btilder * BdotvR))/(Sr - Sstar);

    } else if (dim == 2) {

        double pTl = pl + 0.5 * (Bxl*Bxl + Btildel*Btildel + Bzl*Bzl);
        double pTr = pr + 0.5 * (Bxr*Bxr + Btilder*Btilder + Bzr*Bzr);

        double BdotvL = Bxl*uxl + Btildel*uyl + Bzl*uzl;
        double BdotvR = Bxl*uxr + Btilder*uyr + Bzr*uzr;

        Sstar = (pTr - pTl - Btildel*Btildel + Btilder*Btilder + rhol * uyl * (Sl - uyl) - rhor * uyr * (Sr - uyr)) / (rhol * (Sl - uyl) - rhor * (Sr - uyr)); // eq (10.37)

        pTLstar = rhol*(Sl - uyl)*(Sstar - uyl) + pTl - Btildel*Btildel + Wstarl[6]*Wstarl[6];
        pTRstar = rhor*(Sr - uyr)*(Sstar - uyr) + pTr - Btilder*Btilder + Wstarr[6]*Wstarr[6];

        rhoLstar = rhol * (Sl - uyl)/(Sl - Sstar);
        rhoRstar = rhor * (Sr - uyr)/(Sr - Sstar);

        rhouxLstar = rhol*uxl * (Sl - uyl)/(Sl - Sstar) - (Wstarl[6]*Wstarl[5] - Btildel*Bxl)/(Sl - Sstar);
        rhouxRstar = rhor*uxr * (Sr - uyr)/(Sr - Sstar) - (Wstarr[6]*Wstarr[5] - Btilder*Bxr)/(Sr - Sstar);

        rhouyLstar = rhoLstar * Sstar;
        rhouyRstar = rhoRstar * Sstar;

        rhouzLstar = rhol*uzl * (Sl - uyl)/(Sl - Sstar) - (Wstarl[6]*Wstarl[7] - Btildel*Bzl)/(Sl - Sstar);
        rhouzRstar = rhor*uzr * (Sr - uyr)/(Sr - Sstar) - (Wstarr[6]*Wstarr[7] - Btilder*Bzr)/(Sr - Sstar);

        ELstar = Ql[4] * (Sl - uyl)/(Sl - Sstar) + (pTLstar*Sstar - pTl*uyl - (Wstarl[6] * BhllvhllL - Btildel * BdotvL))/(Sl - Sstar);
        ERstar = Qr[4] * (Sr - uyr)/(Sr - Sstar) + (pTRstar*Sstar - pTr*uyr - (Wstarr[6] * BhllvhllR - Btilder * BdotvR))/(Sr - Sstar);

    }

    Qstarl = {rhoLstar, rhouxLstar, rhouyLstar, rhouzLstar, ELstar, Wstarl[5], Wstarl[6], Wstarl[7]};
    Qstarr = {rhoRstar, rhouxRstar, rhouyRstar, rhouzRstar, ERstar, Wstarr[5], Wstarl[6], Wstarr[7]};

    if (isDivClean) {
        Qstarl[8] = Qlnew[8];
        Qstarr[8] = Qrnew[8];
    }

    for (int n = 0; n < nVar; n++) {
        Fstarl[n] = Fl[n] + Sl * (Qstarl[n] - Qlnew[n]);
        Fstarr[n] = Fr[n] + Sr * (Qstarr[n] - Qrnew[n]);
    }

    if (Sl <=0 && Sstar >= 0) {
        Fhllc = Fstarl;
        return;
    }

    if (Sstar <=0 && Sr >= 0) {
        Fhllc = Fstarr;
        return;
    }
}

/**
 * @brief Force flux
 *
 *        Toro's book pg 600 section 18.2
 *
 * @param   Ql          Left vector of conservative variables
 * @param   Qr          Right vector of conservative variables
 * @param   Fforce      FORCE intercell flux
 * @param   dx          Spatial discretisation length in x-direction
 * @param   dt          Time step
 */
void NumericalFlux::FORCEFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fforce, const double& dx, const double& dt, int dim, const double& ch) {

    ADouble Fl, Fr, Flf, Qtemp, Flw;

    Fl = convertVars.MhdFlux((ADouble) Ql, dim, ch);
    Fr = convertVars.MhdFlux((ADouble) Qr, dim, ch);

    for (int n = 0; n < nVar; n++) {
        Flf[n] = 0.5 * (Fl[n] + Fr[n]) - 0.5 * (dx/dt) * (Qr[n] - Ql[n]);   // Lax–Friedrichs flux, eq (18.6)
        Qtemp[n] = 0.5 * (Ql[n] + Qr[n]) - 0.5 * (dt/dx) * (Fr[n] - Fl[n]); // Two–step Lax–Wendroff flux, eq (18.7)
    }

    Flw = convertVars.MhdFlux(Qtemp, dim, ch); // Two-step Lax–Wendroff = Richtmyer flux, eq (18.7)

    for (int n = 0; n < nVar; n++) {
        Fforce[n] = 0.5 * (Flf[n] + Flw[n]); // FORCE flux, eq (18.8)
    }
}




/** Pressure–Based Wave Speed Estimates (ideal gases) **/
// double ql, qr;
//
// /** Two-rarefaction Riemann solver TRRS for computing Pstar
//  *  Toro's book page 301 & 330, section 9.4.1
//  */
// double z = (Gamma - 1)/(2.0 * Gamma); // eq (10.64)
// double pLR = pow(pl/pr, z);
// double ustar = (pLR*ul/al+ur/ar+2.0*(pLR-1.0)/(Gamma-1.0))/(pLR/al+1.0/ar); // eq (9.35)
// double pstar = 0.5*(pl*pow(1.0+(Gamma-1.0)/(2.0*al)*(ul-ustar),1.0/z)+pr*pow(1.0+(Gamma-1.0)/(2.0*ar)*(ustar-ur),1.0/z)); // eq (9.36)
//
// /** Toro's book, page 330, eq (10.59 & 10.60) **/
// if (pstar <= pl) {
// 	ql = 1.0;
// } else { ql = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pstar/pl-1.0)); }
//
// if (pstar <= pr) {
// 	qr = 1.0;
// } else { qr = sqrt(1.0+(Gamma+1.0)/(2.0*Gamma)*(pstar/pr-1.0)); }
//
// Sl = ul - al * ql;
// Sr = ur + ar * qr;

// Add here different Wave speed estimates ...
/// (10.47)
// Sl = ul - al;
// Sr = ur + ar;

/// (10.48)
// Sl = min(ul - al, ur - ar);
// Sr = max(ul + al, ur + ar);

/// Roe average eigenvalues
// double utilde = (sqrt(rhol) * ul + sqrt(rhor) * ur)/(sqrt(rhol) + sqrt(rhor)); // (10.50)
// double el = computeInternalEnergyFromEoS(rhol, pl);
// double er = computeInternalEnergyFromEoS(rhor, pr);
// double hl = el + pl/rhol;
// double hr = er + pr/rhor;
// double Hl = 0.5 * ul * ul + hl;
// double Hr = 0.5 * ur * ur + hr;
// double Htilde = (sqrt(rhol) * Hl + sqrt(rhor) * Hr)/(sqrt(rhol) + sqrt(rhor)); // (10.51)
// double atilde = sqrt((Gamma - 1) * (Htilde - 0.5 * utilde * utilde));
// Sl = utilde - atilde;
// Sr = utilde + atilde;

/// Einfeldt
// double eta2 = 0.5 * ((sqrt(rhol) * sqrt(rhor))/((sqrt(rhol) + sqrt(rhor)) * (sqrt(rhol) + sqrt(rhor))));
// double dbar = sqrt((sqrt(rhol) * al * al + sqrt(rhor) * ar * ar)/(sqrt(rhol) + sqrt(rhor)) + eta2 * (ur - ul) * (ur - ul));
// double utilde = (sqrt(rhol) * ul + sqrt(rhor) * ur)/(sqrt(rhol) + sqrt(rhor)); // (10.50)
// Sl = utilde - dbar;
// Sr = utilde + dbar;

/// (10.56)
// double Splus = max(fabs(ul) + al, fabs(ur) + ar);
// Sl = - Splus;
// Sr = Splus;