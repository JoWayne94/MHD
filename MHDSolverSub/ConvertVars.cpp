/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : ConvertVars.cpp                            */
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
/*   Description : Convert primitive variables to conservative        */
/*                 variables and vice versa. Convert conservative     */
/*                 variables to MHD flux.                             */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "ConvertVars.h"

using namespace std;

/**
 * @brief Constructor of the ConvertVars class
 */
ConvertVars::ConvertVars() {
}

/**
 * @brief Destructor of the ConvertVars class
 */
ConvertVars::~ConvertVars() {
}

/**
 * @brief Set variables
 *
 * @param   nvar     Number of conservative variables
 * @param   gamma    (Constant) adiabatic index (ratio of specific heats)
 */
void ConvertVars::SetVariables(unsigned int nvar, double gamma, bool isdivclean) {

    nVar = nvar;
    Gamma = gamma;
    isDivClean = isdivclean;

    equationofStates.SetVariables(Gamma);

}

/**
 * @brief  Conservative to primitive routine
 *
 * @return  Primitive variables vector
 */
ADouble ConvertVars::conservativeToprimitive(const ADouble& U) {

    double rho = U[0];
    double u_x = U[1]/U[0];
    double u_y = U[2]/U[0];
    double u_z = U[3]/U[0];
    double E   = U[4];
    double B_x = U[5];
    double B_y = U[6];
    double B_z = U[7];

    double kin = (double) 0.5*rho*(u_x*u_x + u_y*u_y + u_z*u_z);
    double mp  = (double) 0.5*(B_x*B_x + B_y*B_y + B_z*B_z);
    double e   = (double) (E - kin - mp)/rho;
    double p   = equationofStates.computePressureFromEoS((double) rho, (double) e);

    w[0] = rho;
    w[1] = u_x;
    w[2] = u_y;
    w[3] = u_z;
    w[4] = p;
    w[5] = B_x;
    w[6] = B_y;
    w[7] = B_z;
    if (isDivClean) w[8] = U[8]; // psi

    return w;
}

/**
 * @brief  Primitive to conservative routine
 *
 * @return  Conservative variables vector
 */
ADouble ConvertVars::primitiveToconservative(const ADouble& W) {

    double rho = W[0];
    double u_x = W[1];
    double u_y = W[2];
    double u_z = W[3];
    double p   = W[4];
    double B_x = W[5];
    double B_y = W[6];
    double B_z = W[7];

    double kin = (double) 0.5*rho*(u_x*u_x + u_y*u_y + u_z*u_z);
    double mp  = (double) 0.5*(B_x*B_x + B_y*B_y + B_z*B_z);
    double e   = equationofStates.computeInternalEnergyFromEoS((double) rho, (double) p);
    double E   = (double) rho*e+kin+mp;

    Q[0] = rho;
    Q[1] = rho*u_x;
    Q[2] = rho*u_y;
    Q[3] = rho*u_z;
    Q[4] = E;
    Q[5] = B_x;
    Q[6] = B_y;
    Q[7] = B_z;
    if (isDivClean) Q[8] = W[8]; // psi

    return Q;
}

/**
 * @brief  Conservative variables to flux vector
 *
 * @return  MHD flux vector
 */
ADouble ConvertVars::MhdFlux(const ADouble& U, int dim, const double& ch) {

    Wtemp = conservativeToprimitive(U);
    double rho = Wtemp[0];
    double u_x = Wtemp[1];
    double u_y = Wtemp[2];
    double u_z = Wtemp[3];
    double p   = Wtemp[4];
    double B_x = Wtemp[5];
    double B_y = Wtemp[6];
    double B_z = Wtemp[7];
    double E   = U[4];

    double mp  = (double) 0.5*(B_x*B_x + B_y*B_y + B_z*B_z);

    if (dim == 1) {

        F[0] = rho*u_x;
        F[1] = rho*u_x*u_x+p+mp-B_x*B_x;
        F[2] = rho*u_x*u_y     -B_x*B_y;
        F[3] = rho*u_x*u_z     -B_x*B_z;
        F[4] = (E+p+mp)*u_x-(u_x*B_x + u_y*B_y + u_z*B_z)*B_x;
        if (isDivClean) F[5] = Wtemp[8]; // psi
        else F[5] = 0.0;
        F[6] = B_y*u_x - B_x*u_y;
        F[7] = B_z*u_x - B_x*u_z;
        if (isDivClean) F[8] = ch*ch * B_x;

    } else if (dim == 2) {

        F[0] = rho*u_y;
        F[1] = rho*u_y*u_x     -B_y*B_x;
        F[2] = rho*u_y*u_y+p+mp-B_y*B_y;
        F[3] = rho*u_y*u_z     -B_y*B_z;
        F[4] = (E+p+mp)*u_y-(u_x*B_x + u_y*B_y + u_z*B_z)*B_y;
        F[5] = B_x*u_y - B_y*u_x;
        if (isDivClean) F[6] = Wtemp[8];
        else F[6] = 0.0;
        F[7] = B_z*u_y - B_y*u_z;
        if (isDivClean) F[8] = ch*ch * B_y;

    }

    return F;
}