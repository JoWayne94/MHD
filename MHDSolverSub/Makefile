# Declare compiler, required flags, and dependencies
CC = g++
CXXFLAGS = -std=c++17 -pedantic -O2 -g
LDLIBS = -lboost_program_options
HDRS = MHD.h FluxUpdate.h NumericalFlux.h ConvertVars.h EquationofStates.h ExactRPSolver.h TypeDefs.h
TARGET = Solve
# SRCS = MHDSolver.cpp MHD.cpp FluxUpdate.cpp NumericalFlux.cpp ConvertVars.cpp EquationofStates.cpp
OBJS = MHDSolver.o MHD.o FluxUpdate.o NumericalFlux.o ConvertVars.o EquationofStates.o ExactRPSolver.o TypeDefs.o

# Define make target
default: $(TARGET)
all: $(TARGET)

# Generate object files
%.cpp:

%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

# Link object file and produce executable
$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: testBrioWu1D testToro1D testToro2Dx testToro2Dy testToro2Dxy testCylindricalExp testBrioWu2Dx testBrioWu2Dy testBrioWu2Dxy testOrszagTangNoD testOrszagTang testKelvinHelmNoD testKelvinHelm clean cleanPlots cleanData

# Test cases for the code
# N = 100, 200, 400, 800 for exact   0.01, 0.005, 0.0025, 0.00125
testBrioWu1D:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=800 --yBottomDomain=0 --yTopDomain=1 --Nx=800 --Ny=1 --nVar=8 --Gamma=2 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=1 --name="Data/BrioWu1dN800.dat"

testToro1D:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=0.0025 --Nx=400 --Ny=1 --nVar=8 --Gamma=1.4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=2 --name="Data/Sod1dN400.dat"

# 2D, one res
testToro2Dx:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=400 --Ny=400 --nVar=8 --Gamma=1.4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=2 --name="Data/Sod2dxN400.dat"

testToro2Dy:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=400 --Ny=400 --nVar=8 --Gamma=1.4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=3 --name="Data/Sod2dyN400.dat"

testToro2Dxy:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=400 --Ny=400 --nVar=8 --Gamma=1.4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=4 --name="Data/Sod2dxyN400.dat"

testCylindricalExp:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=2 --yBottomDomain=0 --yTopDomain=2 --Nx=101 --Ny=101 --nVar=8 --Gamma=1.4 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=5 --name="Data/CylindricalExp.dat"

testBrioWu2Dx:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=800 --yBottomDomain=0 --yTopDomain=800 --Nx=800 --Ny=800 --nVar=8 --Gamma=2 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=6 --name="Data/BrioWu2dxN800.dat"

testBrioWu2Dy:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=800 --yBottomDomain=0 --yTopDomain=800 --Nx=800 --Ny=800 --nVar=8 --Gamma=2 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=7 --name="Data/BrioWu2dyN800.dat"

testBrioWu2Dxy:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=800 --yBottomDomain=0 --yTopDomain=800 --Nx=800 --Ny=800 --nVar=8 --Gamma=2 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=8 --name="Data/BrioWu2dxyN800.dat"

testOrszagTangNoD:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=256 --Ny=256 --nVar=8 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=9 --name="Data/OTNoD/OrszagTang"

testOrszagTang:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=1024 --Ny=1024 --nVar=9 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=9 --name="Data/OTRef/OrszagTang"

testKelvinHelmNoD:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=-1 --yTopDomain=1 --Nx=256 --Ny=512 --nVar=8 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=10 --name="Data/KelvinHelmNoD/KelvinHelm"

testKelvinHelm:$(TARGET)
	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=-1 --yTopDomain=1 --Nx=256 --Ny=512 --nVar=9 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=10 --name="Data/KelvinHelm/KelvinHelm"

# Clean commands
clean:
	rm -f $(TARGET) *.o

cleanPlots:
	rm -f Plots/*

cleanData:
	rm -f Data/*

#testDaiWoodward:$(TARGET)
#	./$(TARGET) --xLeftDomain=0 --xRightDomain=1 --yBottomDomain=0 --yTopDomain=1 --Nx=200 --Ny=1 --nVar=8 --Gamma=1.6667 --Cfl=0.9 --NxGhost=2 --NyGhost=2 --testNumber=3