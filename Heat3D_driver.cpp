
#include <iostream>

#include "nrutil.h"
#include "util.h"

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
    double*** ic = allocate3DMatrix<double>(1,nx,1,ny,1,nz);
    double* xx = dvector(1,nx);
    double* yy = dvector(1,ny);
    double* zz = dvector(1,nz);

    // create linspace
    createLinespace(xx, 0, l_x, nx);
    createLinespace(yy, 0, l_y, ny);
    createLinespace(zz, 0, l_z, nz);

    createGaussian(0.5, ic, xx, nx, yy, ny, zz, nz);

    // create a boundary condition
    double** xy_bc = dmatrix(1,nx,1,ny);
    double** xz_bc = dmatrix(1,nx,1,nz);
    double** yz_bc = dmatrix(1,ny,1,nz);

    createBC(0,xy_bc,nx,ny);
    createBC(0,xz_bc,nx,nz);
    createBC(0,yz_bc,ny,nz);

    // initialize the solver
    HeatFTCD3D he(nx, ny, nz, l_x, l_y, l_z,alpha, dx, dy, dz, dt);
    std::cout << "Setting IC with nx= " << nx << " ny=" << ny << " nz=" << nz;
    std::cout << std::endl;

    std::cout << "IC[1][1][1]= " << ic[1][1][1] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[1][2][1] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[3][2][1] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[3][3][3] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[4][5][6] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[13][13][13] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[50][50][50] << std::endl;
    std::cout << "IC[nx][ny][nz]= " << ic[100][100][100] << std::endl;

    he.setIC(ic, nx, ny, nz);
    he.setBC(xy_bc, xz_bc, yz_bc, xy_bc, xz_bc, yz_bc, nx, ny, nz);

    // solve
    he.solve(10);

}

