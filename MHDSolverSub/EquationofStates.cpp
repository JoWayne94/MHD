/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : EquationofStates.cpp                       */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 23/01/2022                                                */
/*                                                                    */
/*   Programmer :                                        */
/*                                                                    */
/*   Description : Compute parameters using defined EoS.              */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "EquationofStates.h"
#include <cmath>

using namespace std;

/**
 * @brief Constructor of the EquationofStates class
 */
EquationofStates::EquationofStates() {
}

/**
 * @brief Destructor of the EquationofStates class
 */
EquationofStates::~EquationofStates() {
}

/**
 * @brief Set variables
 *
 * @param   gamma   (Constant) adiabatic index (ratio of specific heats)
 */
void EquationofStates::SetVariables(double gamma) {
    Gamma = gamma;
}

double EquationofStates::computePressureFromEoS(const double& rho, const double& e) const {

    // Ideal Gas EoS
    double p = rho*e*(Gamma - 1.0);

    // Add here other EoS ...

    return p;
}

double EquationofStates::computeInternalEnergyFromEoS(const double& rho, const double& p) const {

    // Ideal Gas EoS
    double e = p/(rho*(Gamma - 1.0));

    // Add here other EoS ...

    return e;
}

double EquationofStates::computeSoundSpeedFromEoS(const double& rho, const double& p) const {

    // Ideal Gas EoS
    double a = sqrt(Gamma*p/rho);

    // Add here other EoS ...

    return a;
}