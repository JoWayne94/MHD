#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "ConvertVars.h"
#include "EquationofStates.h"
#include "TypeDefs.h"

using namespace std;

class NumericalFlux {

public:

    NumericalFlux();
    ~NumericalFlux();

    void SetVariables(unsigned int nvar, double gamma, bool isdivclean);
    void WaveEstimates(const ADouble & Ql, const ADouble & Qr, double& Sl, double& Sr, int dim);
    void HLLFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fhll, int dim, const double& ch);
    void HLLCFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fhllc, int dim, const double& ch);
    void FORCEFlux(const ADouble& Ql, const ADouble& Qr, ADouble& Fforce, const double& dx, const double& dt, int dim, const double& ch);

private:

    /// Classes required
    ConvertVars convertVars;
    EquationofStates equationofStates;

    unsigned int nVar;
    double Sl, Sr, Gamma;
    bool isDivClean;

};

#endif