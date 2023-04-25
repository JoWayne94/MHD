#ifndef EXACTRPSOLVER_H
#define EXACTRPSOLVER_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "TypeDefs.h"
#include "EquationofStates.h"
#include "ConvertVars.h"

using namespace std;

class ExactRPSolver {

public:

    ExactRPSolver();
    ~ExactRPSolver();

    void SetVariables(const ADouble& wl, const ADouble& wr, double gamma, unsigned int nvar);
    void ComputeExactSolution(const double& time, const VDouble& x, const double& xd, VectArr& uExact);

private:

    double Gamma;
    unsigned int nVar;

    double G1,G2,G3,G4,G5,G6,G7;

    double rhol, rhor, ul, ur, pl, pr, al, ar;

    EquationofStates equationofStates;
    ConvertVars convertVars;

    double GuessPressure() const;
    double fl_NR(const double& p) const;
    double fr_NR(const double& p) const;
    double f_NR(const double& p);
    double df_NR(const double& p) const;
    ADouble Sample(const double& pstar, const double& ustar, const int& Pattern, const double& csi) const;

};

#endif