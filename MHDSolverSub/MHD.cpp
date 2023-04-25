/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : MHD.cpp                                    */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 10/02/2022                                                */
/*                                                                    */
/*   Programmer :                                         */
/*                                                                    */
/*   Description : Driver for the MHD solver                          */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "MHD.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <utility>

using namespace std;

/**
 * @brief Constructor for the MHD class
 */
MHD::MHD() {
}

/**
 * @brief Destructor for the MHD class
 *
 * 	      Clean up memory
 */
MHD::~MHD() {
    delete u;
    delete uBc;
    delete wExact;
    delete xCells;
    delete yCells;
}

/**
 * @brief Set domain size
 *
 *	      Initialize Lx and Ly, xLeftDomain, xRightDomain, yBottomDomain, and yTopDomain variables
 *
 * @param   xlength        Length of the domain in the x-direction
 * @param   ylength        Length of the domain in the y-direction
 * @param   xLeft          Left domain boundary in the x-direction
 * @param   xRight         Right domain boundary in the x-direction
 * @param   yBottom        Bottom domain boundary in the y-direction
 * @param   yTop           Top domain boundary in the y-direction
 */
void MHD::SetDomainSize(double xlength, double ylength, double xLeft, double xRight, double yBottom, double yTop) {
    Lx = xlength;
    Ly = ylength;
    xLeftDomain = xLeft;
    xRightDomain = xRight;
    yBottomDomain = yBottom;
    yTopDomain = yTop;
}

/**
 * @brief Set number of cells
 *
 *  	  Initialize Nx and Ny arrays, and other variables
 *
 * @param   nx             Number of cells in the x-direction
 * @param   ny             Number of cells in the y-direction
 * @param   nxGhost        Number of fictitious cells in the x-direction
 * @param   nyGhost        Number of fictitious cells in the y-direction
 */
void MHD::SetNumCells(unsigned int nx, unsigned int ny, unsigned int nxGhost, unsigned int nyGhost) {
    Nx = nx;
    Ny = ny;
    NxGhost = nxGhost;
    NyGhost = nyGhost;
    N = Nx * Ny;
}

/**
 * @brief Set Courant–Friedrichs–Lewy number
 *
 * @param   cfl     Courant–Friedrichs–Lewy number
 */
void MHD::SetCourantNumber(double cfl) {
    Cfl = cfl;
}

/**
 * @brief Set mesh size
 *
 * @param   deltax        Mesh size in the x-direction
 * @param   deltay        Mesh size in the y-direction
 */
void MHD::SetMeshSize(double deltax, double deltay) {
    dx = deltax;
    dy = deltay;
}

/**
 * @brief Set number of conservative variables
 *
 * @param   nvar     Number of conservative variables
 */
void MHD::SetnVar(unsigned int nvar) {
    nVar = nvar;
    isDivClean = false;
    if (nVar == 9) isDivClean = true;
}

/**
 * @brief Set adiabatic index
 *
 * @param   gamma    (Constant) adiabatic index (ratio of specific heats)
 */
void MHD::SetGamma(double gamma) {
    Gamma = gamma;
}

/**
 * @brief Set test number
 *
 * @param   testNum   Choice of pre-defined test cases
 */
void MHD::SetTestNum(unsigned int testNum) {
    testNumber = testNum;
}

/**
 * @brief Set output file name
 *
 * @param   Name   Data .dat file name
 */
void MHD::SetFileName(std::string Name) {
    name = std::move(Name);
}

/**
 * @brief Set BC types
 *
 * @param   xbc   BC type in x-direction
 * @param   ybc   BC type in y-direction
 */
//void MHD::SetBCType(BCtype xbc, BCtype ybc) {
//    xBC = xbc;
//    yBC = ybc;
//}

/**
 * @brief Initialise conservative variable vectors
 *
 *	      Initialise u, uBc, uExact, xCells, and yCells vectors to satisfy initial conditions
 */
void MHD::Initialise() {

    u = new VectArr;
    u->SetSize(Nx, Ny);
    uBc = new VectArr;
    uBc->SetSize(Nx + 2 * NxGhost, Ny + 2 * NyGhost);
    wExact = new VectArr;
    wExact->SetSize(Nx, Ny);
    xCells = new VDouble(N); // x coordinates
    yCells = new VDouble(N); // y coordinates

    convertVars.SetVariables((unsigned int) nVar, (double) Gamma, isDivClean);
}

/**
 * @brief Construction of computational domain
 */
void MHD::ConstructFVDomain(VDouble& xCellCentres, VDouble& yCellCentres) const {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            xCellCentres[i + j*Nx] = xLeftDomain + dx * (i + 0.5);
            yCellCentres[i + j*Nx] = yBottomDomain + dy * (j + 0.5);
        }
    }
}

/**
 * @brief Set initial conditions
 *
 *        Toro's tests 1, 2, 3, 4, 5, 6, 7 (Page 334), and other tests
 */
void MHD::SetInitialConditions(const VDouble& xCellCentres, const VDouble& yCellCentres, VectArr& uIC, double& finalt, double& xDiscontinuity, double& yDiscontinuity) {

    ADouble wl, wr, Ql, Qr;

    switch(testNumber) {
        case(1): // Brio & Wu
        {
            finalt = 80.0;
            xDiscontinuity = 400.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0};

            break;
        }

        case 2 ... 4: // Sod test
        {
            finalt = 0.25;
            xDiscontinuity = 0.5;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

            break;
        }

        case(5): // Cylindrical Explosion
        {
            finalt = 0.25;
            xDiscontinuity = 0.4;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};

            break;
        }

        case(6): // Brio & Wu 2Dx
        {
            finalt = 80.0;
            xDiscontinuity = 400.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0};

            break;
        }

        case(7): // Brio & Wu 2Dy
        {
            finalt = 80.0;
            xDiscontinuity = 400.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, -1.0, 0.75, 0.0};

            break;
        }

        case(8): // Brio & Wu 2Dxy
        {
            finalt = 80.0;
            xDiscontinuity = 400.0;

            xBC = Neumann;
            yBC = Neumann;

            wl = {1.0, 0.0, 0.0, 0.0, 1.0, (-0.75 - 1.0) * cos(45*M_PI/180), (0.75 + 1.0) * cos(45*M_PI/180), 0.0};
            wr = {0.125, 0.0, 0.0, 0.0, 0.1, (-0.75 + 1.0) * cos(45*M_PI/180), (0.75 - 1.0) * cos(45*M_PI/180), 0.0};

            break;
        }

        case(9): // Orszag-Tang with and without div cleaning
        {
            finalt = 1.0;

            xBC = Periodic;
            yBC = Periodic;

            break;
        }

        case(10): // Kelvin-Helmholtz Instability
        {
            finalt = 20.0;

            xBC = Periodic;
            yBC = Reflective;

            break;
        }

        default:
        {
            std::cerr << "Test case number invalid.";
            std::exit(EXIT_FAILURE);
        }

//        case(): // Dai and Woodward
//        {
//            finalt = 0.08;
//            xDiscontinuity = 0.5;
//
//            wl = {1.0, 10.0, 0.0, 0.0, 20.0, 5.0/(2*sqrt(M_PI)), 5.0/(2*sqrt(M_PI)), 0.0};
//            wr = {1.0, -10.0, 0.0, 0.0, 1.0, 5.0/(2*sqrt(M_PI)), 5.0/(2*sqrt(M_PI)), 0.0};
//
//            break;
//        }
    }

    /// Assign initial condition values
    if (testNumber == 3 || testNumber == 7) { // Toro 2Dy, BrioWu 2Dy

        Ql = convertVars.primitiveToconservative((ADouble) wl);
        Qr = convertVars.primitiveToconservative((ADouble) wr);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( yCellCentres[i + j*Nx] < xDiscontinuity ) { uIC(i, j) = Ql; }
                else { uIC(i, j) = Qr; }

            }
        }

    } else if (testNumber == 4 || testNumber == 8) { // Toro 2Dxy, BrioWu 2Dxy

        Ql = convertVars.primitiveToconservative((ADouble) wl);
        Qr = convertVars.primitiveToconservative((ADouble) wr);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] + yCellCentres[i + j*Nx] < 800.0 ) { uIC(i, j) = Ql; } // 1.0
                else { uIC(i, j) = Qr; }

            }
        }

    } else if (testNumber == 5) { // Cylindrial explosion

        Ql = convertVars.primitiveToconservative((ADouble) wl);
        Qr = convertVars.primitiveToconservative((ADouble) wr);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( sqrt((xCellCentres[i + j*Nx] - 1)*(xCellCentres[i + j*Nx] - 1) + (yCellCentres[i + j*Nx] - 1)*(yCellCentres[i + j*Nx] - 1)) < xDiscontinuity ) { uIC(i, j) = Ql; }
                else { uIC(i, j) = Qr; }

            }
        }

    } else if (testNumber == 9) {

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                wl[0] = ((double) 5/3)*((double) 5/3);
                wl[1] = -sin(2.0 * M_PI * yCellCentres[i + j*Nx]);
                wl[2] = sin(2.0 * M_PI * xCellCentres[i + j*Nx]);
                wl[3] = 0.0;
                wl[4] = ((double) 5/3);
                wl[5] = -sin(2.0 * M_PI * yCellCentres[i + j*Nx]);
                wl[6] = sin(4.0 * M_PI * xCellCentres[i + j*Nx]);
                wl[7] = 0.0;
                if (isDivClean) wl[8] = 0.0;

                uIC(i, j) = convertVars.primitiveToconservative((ADouble) wl);
            }
        }

    } else if (testNumber == 10) {

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                wl[0] = 1.0;
                wl[1] = 0.5 * tanh(20.0 * yCellCentres[i + j*Nx]);
                wl[2] = 0.01 * sin(2.0 * M_PI * xCellCentres[i + j*Nx]) * exp(- yCellCentres[i + j*Nx]*yCellCentres[i + j*Nx] / (0.1 * 0.1));
                wl[3] = 0.0;
                wl[4] = ((double) 3/5);
                wl[5] = 0.1 * sqrt(1.0) * cos((double) M_PI / 3.0);
                wl[6] = 0.0;
                wl[7] = 0.1 * sqrt(1.0) * sin((double) M_PI / 3.0);
                if (isDivClean) wl[8] = 0.0;

                uIC(i, j) = convertVars.primitiveToconservative((ADouble) wl);
            }
        }

    } else { // BrioWu 1D, Toro 1D, Toro 2Dx, BrioWu 2Dx

        Ql = convertVars.primitiveToconservative((ADouble) wl);
        Qr = convertVars.primitiveToconservative((ADouble) wr);

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if ( xCellCentres[i + j*Nx] < xDiscontinuity ) { uIC(i, j) = Ql; }
                else { uIC(i, j) = Qr; }

            }
        }

    }

    /// Save left and right states of the RP
    RP_LeftState = wl;
    RP_RightState = wr;
}

/**
 * @brief Set boundary conditions
 */
void MHD::SetBoundaryConditions(const VectArr& uIC, VectArr& ubc) const {

    /// Copy the valid state
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            ubc(i + NxGhost, j + NyGhost) = uIC(i, j);
        }
    }

    switch(xBC)
    {
        case(Neumann):
        {
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < NxGhost; i++) {
                    ubc(i, j + NyGhost) = uIC(0, j);
                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(Nx - 1, j);
                }
            }
            break;
        }
        case(Periodic):
        {
            /// Periodic BC in x-direction
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < NxGhost; i++) {
                    ubc(i, j + NyGhost) = uIC(Nx - NxGhost + i, j);
                    ubc(Nx + NxGhost + i, j + NyGhost) = uIC(NxGhost + i, j);
                }
            }
            break;
        }
        default:
        {
            std::cerr << "Boundary conditions not recognized.";
            std::exit(EXIT_FAILURE);
        }
    }

    switch(yBC)
    {
        case(Neumann):
        {
            /// Neumann boundary conditions
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {
                    ubc(i + NxGhost, j) = uIC(i, 0);
                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, Ny - 1);
                }
            }
            break;
        }
        case(Periodic):
        {
            /// Implement here other boundary conditions ... (reflective for example)
            /// Periodic BC in y-direction
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {
                    ubc(i + NxGhost, j) = uIC(i, Ny - NyGhost + j);
                    ubc(i + NxGhost, Ny + NyGhost + j) = uIC(i, NyGhost + j);
                }
            }
            break;
        }
        case(Reflective):
        {
            /// Reflective BC in y-direction
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < NyGhost; j++) {

                    for (int k = 0; k < nVar; k++) {
                        /// Reflecting u_y and B_y
                        if (k == 2 || k == 6) {
                            ubc(i + NxGhost, j, k) = -uIC(i, NyGhost - 1 - j, k);
                            ubc(i + NxGhost, Ny + NyGhost + j, k) = -uIC(i, Ny - 1 - j, k);
                        } else {
                            ubc(i + NxGhost, j, k) = uIC(i, NyGhost - 1 - j, k);
                            ubc(i + NxGhost, Ny + NyGhost + j, k) = uIC(i, Ny - 1 - j, k);
                        }
                    }

                }
            }
            break;
        }
        default:
        {
            std::cerr << "Boundary conditions not recognized.";
            std::exit(EXIT_FAILURE);
        }
    }
}

/**
 * @brief Solution to the user-defined problem to compute conservative variable values at time t + dt
 *        Update all the steps together using given EoS to solve for the conservative variables iteratively
 * 	      Create a new RiemannSolver instance and solve the Riemann problem
 */
void MHD::Update() {

    /// Construct computational domain
    cout << "Constructing spatial domain... " << endl;
    ConstructFVDomain((*xCells), (*yCells));

    /// Set initial conditions
    cout << "Setting initial conditions... " << endl;
    SetInitialConditions((*xCells), (*yCells), (*u), T, RP_xDisc, RP_yDisc);

    /// Initialise time loop
    int iter = 0;
    time = 0.0;
    double ch; // max wave speed

    /// Initialise divergence error and timing output dat files
//    std::remove("Data/KelvinHelm/L1diverror.dat"); // delete file
//    std::ofstream outErr("Data/KelvinHelm/L1diverror.dat");
//
//    std::remove("Data/KelvinHelm/Timing.dat"); // delete file
//    std::ofstream outTime("Data/KelvinHelm/Timing.dat");

    /// FluxUpdate routines
    fluxUpdate.SetVariables((unsigned int) Nx, (unsigned int) Ny, (double) Cfl, (unsigned int) NxGhost, (unsigned int) NyGhost,
                             (unsigned int) nVar, (double) T, (double) Gamma, isDivClean);
    fluxUpdate.Initialise();

    /// Timing parameters
    clock_t start, end;
    double elapsed_time = 0.0;
    double errorTime = 0.0;
    double outputTime = 0.0;
    PairDouble timerxUpdate, timeryUpdate; // dimensional split timings

    cout << "Starting time loop... " << endl;
    do {
        /// Dimensional split x-direction
        // Apply Boundary Conditions
        SetBoundaryConditions((*u), (*uBc));

        /// Compute time step
        ch = 0.0;  // Re-initialise max wave speed for every iteration
        start = clock();
        fluxUpdate.ComputeDt((*uBc), dx, dy, time, dt, ch);
        end = clock();
        elapsed_time = (end - start)/(double)CLOCKS_PER_SEC ;

        // double divB;
        // double L1err = 0.0;
        // double maxerr = 0.0;

//        if (outputTime - dt < time && time < outputTime + dt) {
//
//            GenerateData((*xCells), (*yCells), (*u), (*wExact));
//
//            outputTime += 1.0;
//        }

//        if (iter % 20 == 0) {
//            std::ofstream outIter;
//
//            outIter.open(name + std::to_string(iter) + ".dat");
//
//            outIter << "# This file is the output of MHD.cpp" << std::endl;
//
//            outIter << "# " << std::endl;
//
//            VectArr prim;
//            prim.SetSize(Nx, Ny); // primitive for output
//
//            equationofStates.SetVariables(Gamma);
//            double psi;
//
//            for (int i = 0; i < Nx; i++) {
//                for (int j = 0; j < Ny; j++) {
//                    prim(i, j) = convertVars.conservativeToprimitive((*u)(i, j));
//
//                    outIter << (*xCells)[i + j*Nx] << " " << (*yCells)[i + j*Nx]
//                            << " " << prim(i, j, 0) << " " << prim(i, j, 1) << " " << prim(i, j, 2) << " " << prim(i, j, 3) << " " << equationofStates.computeInternalEnergyFromEoS(prim(i, j, 0), prim(i, j, 4))
//                            << " " << prim(i, j, 4) << " " << prim(i, j, 5) << " " << prim(i, j, 6) << " " << prim(i, j, 7)
//                            << " " << prim(i, j, 8) << " " << sqrt(prim(i, j, 5)*prim(i, j, 5) + prim(i, j, 6)*prim(i, j, 6))/prim(i, j, 7) << "\n";
//                }
//                outIter << "\n";
//            }
//
//            outIter.close();
//        }
//
//        if (errorTime - dt < time && time < errorTime + dt) {
//            for (int i = NxGhost; i < Nx + NxGhost; i++) {
//                for (int j = NyGhost; j < Ny + NyGhost; j++) {
//                    divB = ((*uBc)(i + 1, j, 5) - (*uBc)(i - 1, j, 5))/(2 * dx) + ((*uBc)(i, j + 1, 6) - (*uBc)(i, j - 1, 6))/(2 * dy);
//                    L1err += fabs(divB);
//                    // maxerr = std::max(maxerr, divB);
//                }
//            }
//
//            outErr << time << " " << L1err/N << std::endl;
//            errorTime += 1.0;
//        }

        // Compute numerical fluxes and update solution
        fluxUpdate.UpdatewithFluxes((*uBc), (*u), dx, dt, 1, ch, timerxUpdate);

        /// Dimensional split y-direction
        // Apply Boundary Conditions
        SetBoundaryConditions((*u), (*uBc));

        // Compute numerical fluxes and update solution
        fluxUpdate.UpdatewithFluxes((*uBc), (*u), dy, dt, 2, ch, timeryUpdate);

        /// Operator splitting to handle source term for mixed divergence cleaning
        if (isDivClean) {
            SetBoundaryConditions((*u), (*uBc));

            fluxUpdate.OperatorSplitting((*uBc), (*u), dt, ch);
        }

        iter += 1;
        time += dt;
        std::cout << "Iteration # " << iter << " Time = " << time << std::endl;

        // outTime << iter << " " << elapsed_time << " " << timerxUpdate.first << " " << timerxUpdate.second << " " << timeryUpdate.first << " " << timeryUpdate.second << std::endl;

    } while (time < T);

//    elapsed_time /= iter;
//    cout << "Average time taken for compute dt function in ms is: " << elapsed_time << endl;

    /// Compute Exact Solution
    if (testNumber == 2) {
        exactRPSolver.SetVariables(RP_LeftState, RP_RightState, Gamma, nVar);
        exactRPSolver.ComputeExactSolution(time, (*xCells), RP_xDisc, (*wExact));
    }

    /// Write to output file
    cout << "Writing data file... " << endl;
    GenerateData((*xCells), (*yCells), (*u), (*wExact));
}

/**
 *  @brief  Output x & y coordinates and conservative variable values to dat file
 */
void MHD::GenerateData(const VDouble& xCellCentres, const VDouble& yCellCentres, const VectArr& uOut, const VectArr& wExactOut) {

    std::ofstream outFile;

    outFile.open(name); // + "T" + std::to_string(time) + ".dat"

    outFile << "# This file is the output of MHD.cpp" << std::endl;

    outFile << "# " << std::endl;

    VectArr prim;
    double psi, rhoExact, uExact, pExact, eExact;
    prim.SetSize(Nx, Ny); // primitive for output

    equationofStates.SetVariables(Gamma);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            prim(i, j) = convertVars.conservativeToprimitive(uOut(i, j));

            if (isDivClean) psi = prim(i, j, 8);
            else psi = 0.0;

//            if (testNumber == 2) {
//                rhoExact = wExactOut(i, j, 0);
//                uExact = wExactOut(i, j, 1);
//                pExact = wExactOut(i, j, 4);
//                eExact = equationofStates.computeInternalEnergyFromEoS(rhoExact, pExact);
//            } else {
//                rhoExact = 0.0;
//                uExact = 0.0;
//                pExact = 0.0;
//                eExact = 0.0;
//            }
            if (xCellCentres[i + j*Nx] == yCellCentres[i + j*Nx]) {
                outFile << xCellCentres[i + j*Nx] << " " << yCellCentres[i + j*Nx]
                        << " " << prim(i, j, 0) << " " << prim(i, j, 1) << " " << prim(i, j, 2) << " " << prim(i, j, 3) << " " << equationofStates.computeInternalEnergyFromEoS(prim(i, j, 0), prim(i, j, 4))
                        << " " << prim(i, j, 4) << " " << prim(i, j, 5) << " " << prim(i, j, 6) << " " << prim(i, j, 7)
                        << " " << psi << " " << sqrt(prim(i, j, 5)*prim(i, j, 5) + prim(i, j, 6)*prim(i, j, 6))/prim(i, j, 7) << "\n";
            }
            // << " " << rhoExact << " " << uExact << " " << pExact << " " << eExact << "\n";
        }
        // outFile << "\n";
    }

    outFile.close();
}
