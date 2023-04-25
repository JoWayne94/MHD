#ifndef MHD_H
#define MHD_H
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "ConvertVars.h"
#include "FluxUpdate.h"
#include "TypeDefs.h"
#include "ExactRPSolver.h"
#include "EquationofStates.h"

using namespace std;

class MHD {

public:

    MHD();
    ~MHD();

    void SetDomainSize(double xlength, double ylength, double xLeft, double xRight, double yBottom, double yTop);
    void SetNumCells(unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost);
    void SetCourantNumber(double cfl);
    void SetMeshSize(double deltax, double deltay);
    void SetnVar(unsigned int nvar);
    void SetGamma(double gamma);
    void SetTestNum(unsigned int testNum);
    void SetFileName(std::string Name);
    // void SetBCType(BCtype xbc, BCtype ybc);

    /// Main solver class functions
    void Initialise();
    void Update();

private:

    /// Static member data
    double Cfl, dx, dy, Gamma, time;
    int Nx, Ny, NxGhost, NyGhost, N;
    unsigned int nVar, testNumber;
    double Lx, Ly, xLeftDomain, xRightDomain, yBottomDomain, yTopDomain;
    std::string name;
    BCtype xBC, yBC;
    bool isDivClean;

    VectArr *u = nullptr;
    VectArr *uBc = nullptr;
    VectArr *wExact = nullptr;
    VDouble *xCells = nullptr;
    VDouble *yCells = nullptr;

    ADouble RP_LeftState, RP_RightState;

    /// To be defined parameters
    double dt, T;
    double RP_xDisc, RP_yDisc;

    /// Private functions to set initial and boundary conditions
    void ConstructFVDomain(VDouble& xCellCentres, VDouble& yCellCentres) const ;
    void SetInitialConditions(const VDouble& xCellCentres, const VDouble& yCellCentres, VectArr& uIC, double& finalt, double& xDiscontinuity, double& yDiscontinuity) ;
    void SetBoundaryConditions(const VectArr& uIC, VectArr& ubc) const ;
    void GenerateData(const VDouble& xCellCentres, const VDouble& yCellCentres, const VectArr& uOut, const VectArr& wExactOut) ;

    /// Classes required
    ConvertVars convertVars;
    FluxUpdate fluxUpdate;
    ExactRPSolver exactRPSolver;
    EquationofStates equationofStates;
};

#endif