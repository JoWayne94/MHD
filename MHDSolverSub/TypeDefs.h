#ifndef TYPEDEFS_H
#define TYPEDEFS_H
#pragma once
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>

typedef std::vector<double> VDouble;
typedef std::vector<std::vector<double> > VVDouble;
typedef std::vector<std::vector<std::vector<double> > > VVVDouble;
typedef std::pair<double, double> PairDouble;
typedef std::vector<std::pair<int, int> > VectPair;
typedef std::array<double, 9 > ADouble;

enum BCtype {Neumann=0, Periodic=1, Reflective=2};

using namespace std;

class VectArr{

public:

    VectArr();

    void SetSize(unsigned int nx, unsigned int ny);
    ADouble& operator()(int i, int j);
    const ADouble& operator()(int i, int j) const ;
    double& operator()(int i, int j, int k);
    const double& operator()(int i, int j, int k) const;

private:

    size_t Nx;
    size_t Ny;
    size_t N;

    std::vector<ADouble > m_data; // 2 vectors as one for index and another for data

};

#endif