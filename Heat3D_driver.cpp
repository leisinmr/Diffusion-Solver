
#include <iostream>

#include "Timer.h"
#include "Heat3D.h"



int main(int argc, char** argv)
{
    // vars for the 3D heat equation
    int nx = 100;
    int ny = 100;
    int nz = 100;
    int l_x = 1;
    int l_y = 1;
    int l_z = 1;
    double alpha = 0.001;
    double dx = (double)l_x/nx;
    double dy = (double)l_y/ny;
    double dz = (double)l_z/nz;
    double dt = 0.0005;
   

    // create an initial condition
    Matrix3D ic(nx, ny, nz);
    Vector xx(nx);
    Vector yy(ny);
    Vector zz(nz);

    GaussianFactory g;
    g(ic, 0.5, xx, yy, zz);

    // create a boundary condition


    Heat3D<Matrix3D, double**, Vector> he(l_x, l_y, l_z, nx, ny, nz, alpha, dx, dy, dz, dt);
    
    he.setIC(ic);


}

