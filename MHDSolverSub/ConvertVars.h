#ifndef CONVERTVARS_H
#define CONVERTVARS_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "EquationofStates.h"
#include "TypeDefs.h"

using namespace std;

class ConvertVars {

public:

    ConvertVars();
    ~ConvertVars();

    void SetVariables(unsigned int nvar, double gamma, bool isdivclean);
    ADouble conservativeToprimitive(const ADouble& U);
    ADouble primitiveToconservative(const ADouble& W);
    ADouble MhdFlux(const ADouble& U, int dim, const double& ch);

    bool isDivClean;

private:

    /// Classes required
    EquationofStates equationofStates;

    unsigned int nVar;
    double Gamma;

    ADouble Q;
    ADouble w;
    ADouble F;
    ADouble Wtemp;

};

#endif