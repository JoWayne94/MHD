#ifndef FLUXUPDATE_H
#define FLUXUPDATE_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <initializer_list>
#include "NumericalFlux.h"
#include "ConvertVars.h"
#include "EquationofStates.h"
#include "TypeDefs.h"
#include <ctime>

using namespace std;

class FluxUpdate {

public:

    FluxUpdate();
    ~FluxUpdate();

    void SetVariables(unsigned int nx, unsigned int ny, double cfl, unsigned int nxGhost, unsigned int nyGhost, unsigned int nvar, double finalt, double gamma, bool isdivclean);
    void Initialise();
    void ComputeDt(const VectArr& ubc, const double& deltax, const double& deltay, const double& time, double& deltat, double& Smax);
    void UpdatewithFluxes(const VectArr& ubcOld, VectArr& unew, const double& deltax, const double& deltat, int dim, const double& ch, PairDouble& timerPair);
    void OperatorSplitting(const VectArr& ubcOld, VectArr& unew, const double& deltat, const double& ch) const;

private:

    /// Classes required
    NumericalFlux numericalFlux;
    ConvertVars convertVars;
    EquationofStates equationofStates;

    /// Timing parameters
    clock_t start, end;
    double elapsed_timehts, elapsed_timeflux;

    int Nx, Ny, NxGhost, NyGhost;
    unsigned int nVar;
    double Cfl, T, Gamma;
    VectArr *Ql = nullptr;
    VectArr *Qr = nullptr;

    void HalfTimeStepUpdate(const VectArr& ubc, const double& deltax, const double& deltat, ADouble& qhtsl, ADouble& qhtsr, int k, int l, int dim, const double& ch);

};

#endif