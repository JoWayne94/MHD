#include"TypeDefs.h"

VectArr::VectArr() {
}

void VectArr::SetSize(unsigned int nx, unsigned int ny)
{
    Nx = nx;
    Ny = ny;
    N = Nx * Ny;
    m_data.resize(N);
}

ADouble& VectArr::operator()(int i, int j)
{
    return m_data[i + j*Nx];
}

const ADouble& VectArr::operator()(int i, int j) const
{
    return m_data[i + j*Nx];
}

double& VectArr::operator()(int i, int j, int k)
{
    return m_data[i + j*Nx][k];
}

const double& VectArr::operator()(int i, int j, int k) const
{
    return m_data[i + j*Nx][k];
}
