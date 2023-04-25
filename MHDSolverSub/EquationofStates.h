#ifndef EQUATIONOFSTATES_H
#define EQUATIONOFSTATES_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "TypeDefs.h"

using namespace std;

class EquationofStates {

public:

    EquationofStates();
    ~EquationofStates();

    void SetVariables(double gamma);
    double computePressureFromEoS(const double& rho, const double& e) const;
    double computeInternalEnergyFromEoS(const double& rho, const double& p) const;
    double computeSoundSpeedFromEoS(const double& rho, const double& p) const;

private:

    double Gamma;

};

#endif