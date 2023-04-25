/* -------------------------------------------------------------------*/
/*                                                                    */
/*          Finite Volume Schemes for System of Eqs                   */
/*                                                                    */
/*   Name of the program : MHDSolver.cpp                              */
/*                                                                    */
/*   Purpose : Solve numerically one/two-dimensional PDEs             */
/*             using a first or second order finite volume method.    */
/*             Centred and upwind numerical fluxes are implemented    */
/*             and explored.                                          */
/*                                                                    */
/*   Date : 10/02/2022                                                */
/*                                                                    */
/*   Programmer :                                        */
/*                                                                    */
/*   Description : Controls the input from command line,              */
/*                 storing it in predefined variables. Accept options */
/*                 and sets up the solver for user-defined problem.   */
/*                                                                    */
/* -------------------------------------------------------------------*/

#include "MHD.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <exception>
#include <cmath>
#include <iterator>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char **argv) {

    /// Set up command line options and defaults. Add input options with appropriate default values to run the code
    po::options_description opts("Solves user-defined problem for the given inputs");
    opts.add_options()
            ("help", "Prints help message.")
            ("xLeftDomain",   po::value<double>()      ->default_value(0.0), "Left domain boundary in x-direction. Default = 0.0")
            ("xRightDomain",  po::value<double>()      ->default_value(1.0), "Right domain boundary in x-direction. Default = 1.0")
            ("yBottomDomain", po::value<double>()      ->default_value(0.0), "Bottom domain boundary in y-direction. Default = 0.0")
            ("yTopDomain",    po::value<double>()      ->default_value(1.0), "Top domain boundary in y-direction. Default = 1.0")
            ("Nx",            po::value<unsigned int>()->default_value(100), "Number of cells in x-direction. Default = 100")
            ("Ny",            po::value<unsigned int>()->default_value(100), "Number of cells in y-direction. Default = 100")
            ("Cfl",           po::value<double>()      ->default_value(1),   "Courant–Friedrichs–Lewy number. Default = 1")
            ("nVar",          po::value<unsigned int>()->default_value(4),   "Number of conservative variables. Default = 4")
            ("NxGhost",       po::value<unsigned int>()->default_value(1),   "Number of fictitious cells in x-direction. Default = 1")
            ("NyGhost",       po::value<unsigned int>()->default_value(1),   "Number of fictitious cells in y-direction. Default = 1")
            ("Gamma",         po::value<double>()      ->default_value(((double) 5/3)), "Adiabatic index. Default = 5/3")
            ("testNumber",    po::value<unsigned int>()->default_value(1),   "Test case number. Default = 1")
            ("name",          po::value<std::string>() ->default_value("Data/MHDFV_2dOutput.dat"), "Output file name");
//            ("xBC",           po::value<BCtype>()      ->default_value(Neumann), "Boundary condition type in x-direction. Default is zero Neumann BC")
//            ("yBC",           po::value<BCtype>()      ->default_value(Neumann), "Boundary condition type in y-direction. Default is zero Neumann BC");

    /// Tell boost to parse command-line arguments using list of possible options and generate a map of options and values
    po::variables_map vm;
    /// Parse command-line arguments and store in buffer vm
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    /// Check if user called "--help" option and print usage
    if (vm.count("help")) {
        cout << "Help option called... " << endl;
        cout << opts << endl;
        return 0;
    }

    /// Read variables from vm or command line. Extract options and save as variables
    const double xLeftDomain    = vm["xLeftDomain"].as<double>();
    const double xRightDomain   = vm["xRightDomain"].as<double>();
    const double yBottomDomain  = vm["yBottomDomain"].as<double>();
    const double yTopDomain     = vm["yTopDomain"].as<double>();
    const unsigned int Nx       = vm["Nx"].as<unsigned int>();
    const unsigned int Ny       = vm["Ny"].as<unsigned int>();
    const double Cfl            = vm["Cfl"].as<double>();
    const unsigned int nVar     = vm["nVar"].as<unsigned int>();
    const unsigned int NxGhost  = vm["NxGhost"].as<unsigned int>();
    const unsigned int NyGhost  = vm["NyGhost"].as<unsigned int>();
    const double Gamma          = vm["Gamma"].as<double>();
    const unsigned int testNumber  = vm["testNumber"].as<unsigned int>();
    const std::string name      = vm["name"].as<std::string>();
//    const BCtype xBC            = vm["xBC"].as<BCtype>();
//    const BCtype yBC            = vm["yBC"].as<BCtype>();

    const double Lx = fabs(xLeftDomain - xRightDomain);
    const double Ly = fabs(yBottomDomain - yTopDomain);

    /// Display chosen parameter values
    cout << "Selected length of domain in x-direction = " << Lx << endl;
    cout << "Selected length of domain in y-direction = " << Ly << endl;
    cout << "Selected number of cells in x-direction = " << Nx << endl;
    cout << "Selected number of cells in y-direction = " << Ny << endl;
    cout << "Selected Courant–Friedrichs–Lewy number = " << Cfl << endl;
    cout << "Selected number of conservative variables = " << nVar << endl;
    cout << "Selected number of fictitious cells in x-direction = " << NxGhost << endl;
    cout << "Selected number of fictitious cells in y-direction = " << NyGhost << endl;
    cout << "Selected adiabatic index = " << Gamma << endl;
    cout << "Selected test case = " << testNumber << endl;
//    cout << "Selected boundary condition type in x-direction is " << xBC << endl;
//    cout << "Selected boundary condition type in y-direction is " << yBC << endl;
    cout << endl;

    /// Perform checks on input variables. Check if domain length is valid (positive)
    if ((Lx <= 0) || (Ly < 0)) {
        cout << "Error: Minimum length of the domain must be > 0." << endl;
        return -1;
    }

    /// Check if discretisation is valid (more than or equals to 3 for second order central difference scheme)
//    if ((Nx < 3) || (Ny < 3)) {
//        cout << "Error: Minimum number of cells in x and y-direction is 3. Please pick a value >= 3." << endl;
//        return -1;
//    }

    /// Check if CFL number is valid
    if ((Cfl <= 0) || (Cfl > 1)) {
        cout << "Error: CFL must be 0 < CFL <= 1." << endl;
        return -1;
    }

    /// Check if Gamma is positive
    if (Gamma <= 0) {
        cout << "Error: Minimum value of Gamma must be > 0." << endl;
        return -1;
    }

    /// Check if nVar, NxGhost and NyGhost are positive
    if ((nVar <= 0) || (NxGhost <= 0) || (NyGhost <= 0)) {
        cout << "Error: Minimum values of nVar, NxGhost and/or NyGhost must be > 0." << endl;
        return -1;
    }

    /// Calculate spatial discretisation length in x and y-direction
    const double dx = Lx / ((double) Nx);
    const double dy = Ly / ((double) Ny);

    /// Display dx and dy
    cout << "Mesh size in the x-direction is " << dx << endl;
    cout << "Mesh size in the y-direction is " << dy << endl;
    cout << endl;

    /// Check if chosen number of conservative variables agree with dimensions
    if (Ly > 0 && nVar <= 3) {
        cout << "Error: Number of conservative variables has to be more than 3 for 2D simulation. Please choose a nVar > 3 to meet requirements." << endl;
        return -1;
    }

    /// Creates new instance of the MHD class
    auto* solver = new MHD();

    /// Configure solver
    solver->SetDomainSize(Lx, Ly, xLeftDomain, xRightDomain, yBottomDomain, yTopDomain);
    solver->SetNumCells(Nx, Ny, NxGhost, NyGhost);
    solver->SetCourantNumber(Cfl);
    solver->SetMeshSize(dx, dy);
    solver->SetnVar(nVar);
    solver->SetGamma(Gamma);
    solver->SetTestNum(testNumber);
    solver->SetFileName(name);
    // solver->SetBCType(xBC, yBC);

    /// Run the solver
    if (Ly == 0) {
        cout << "Running 1D simulation... " << endl;
    } else if (Ly > 0) {
        cout << "Running 2D simulation... " << endl;
    }

    solver->Initialise();
    solver->Update();

    return 0;
}