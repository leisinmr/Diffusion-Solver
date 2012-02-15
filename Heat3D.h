
#include <math.h>

#include "util.h"

void createLinespace(double *v, int start, int end, int numIntervals)
{
    std::cout << start << std::endl;
    std::cout << end << std::endl;
    std::cout << numIntervals << std::endl;

    double increment = ((double)end - start) / numIntervals;
    std::cout << increment << std::endl;
    for (int i = 1; i <= numIntervals; ++i)
    {
        v[i] = increment*i;
        std::cout << v[i] << std::endl;
    }
}


void createGaussian(double mu, double ***m, double* x, long xSize, 
                    double* y, long ySize, double *z, long zSize)
{
    for (int i = 1; i <= xSize; ++i)
    {
        for (int j = 1; j <= ySize; ++j)
        {
            for (int k = 1; k <= zSize; ++k)
            {
                m[i][j][k] = exp(-1*pow(5*x[i] - 2.5, 2)) *
                    exp(-1*pow(5*y[i] - 2.5, 2)) *
                    exp(-1*pow(5*z[i] - 2.5, 2));
            }
        }
    }
}

void createBC(int c, double **m, long xSize, long ySize)
{
    for (int i = 1; i <= xSize; ++i)
        for (int j = 1; j <= ySize; ++j)
            m[i][j] = c;
}



/*
 * FTCD 3D Heat Diffusion Solver
 *
 *  
 *
 */
class HeatFTCD3D
{
    typedef double** BC;
    typedef double*** Matrix;
public:

    HeatFTCD3D(int _nx, int _ny, int _nz, int _l_x, int _l_y, int _l_z, 
           double _dx, double _dy, double _dz, double _dt, double alpha) 
        : T(NULL), Tnew(NULL), nx(_nx), ny(_ny), nz(_nz), 
          l_x(_l_x), l_y(_l_y), l_z(_l_z), dx(_dx), dy(_dy), dz(_dz), dt(_dt)
    {
        // allocate space for arrays
        T = allocate3DMatrix<double>(1,nx,1,ny,1,nz);
        Tnew = allocate3DMatrix<double>(1,nx,1,ny,1,nz);

        xy_0 = dmatrix(1,nx,1,ny);
        xz_0 = dmatrix(1,nx,1,nz);
        yz_0 = dmatrix(1,ny,1,nz);
        xy_N = dmatrix(1,nx,1,ny);
        xz_N = dmatrix(1,nx,1,nz);
        yz_N = dmatrix(1,ny,1,nz);

        // solve for constants
        Cx = alpha/pow(dx,2);
        Cy = alpha/pow(dy,2);
        Cz = alpha/pow(dz,2);
    }  

    void setIC(double*** ic, long _x, long _y, long _z)
    {

        std::cout << "Settting IC (in func)..." << std::endl;
        std::cout << "IC[1][1][1]= " << ic[1][1][1] << std::endl;
        for (int i = 1; i <= _x; ++i)
        {
            for (int j = 1; j <= _y; ++j)
            {
                for (int k = 1; k <= _z; ++k)
                {
                    T[i][j][k] = ic[i][j][k];
                }
            }
        }
    }

    void setBC(BC _xy_0, BC _xz_0, BC _yz_0, BC _xy_N, BC _xz_N, BC _yz_N, 
            int xSize, int ySize, int zSize)
    {
        std::cout << "Settting BCs (in func)..." << std::endl;
        // set the XY plane
        for (int i = 1; i <= xSize; ++i)
        {
            for (int j = 1; j <= ySize; ++j)
            {
                xy_0[i][j] = _xy_0[i][j];
                xy_N[i][j] = _xy_N[i][j];
            }
        }
        // set the YZ plane
        for (int j = 1; j <= ySize; ++j)
        {
            for (int k = 1; k <= zSize; ++k)
            {
                yz_0[j][k] = _yz_0[j][k];
                yz_N[j][k] = _yz_N[j][k];
            }
        }
        // set the XZ plane
        for (int i = 1; i <= xSize; ++i)
        {
            for (int k = 1; k <= zSize; ++k)
            {
                xz_0[i][k] = _xz_0[i][k];
                xz_N[i][k] = _xz_N[i][k];
            }
        }
    }

    void solve(int steps)
    {
        /*
         * T_new = T
         */

        // solve the FTCD for a specified # of steps
        for (int n = 0; n < steps; ++n)
        {
            // solve
            solveStep(n);

            // write the output for that step
            
        }
    }

private:
    void resetBC()
    {
        for (int i = 1; i <= nx; ++i)
        {
            for (int j = 1; j <= ny; ++j)
            {
                T[i][j][1]  = xy_0[i][j];
                T[i][j][nz] = xy_N[i][j];
            }
        }
        // set the YZ plane
        for (int j = 1; j <= ny; ++j)
        {
            for (int k = 1; k <= nz; ++k)
            {
                T[1][j][k]  = yz_0[j][k];
                T[nx][j][k] = yz_N[j][k];
            }
        }
        // set the XZ plane
        for (int i = 1; i <= nx; ++i)
        {
            for (int k = 1; k <= nz; ++k)
            {
                T[i][1][k]  = xz_0[i][k];
                T[i][ny][k] = xz_N[i][k];
            }
        }
    }

    void solveStep(int theStep)
    {
        std::cout << "Solving step #"<< theStep << std::endl;
        // temp vars to store the dimensional components
        double xTerm = 0.0, yTerm = 0.0, zTerm = 0.0, tnTerm = 0.0;

        // Go thru each dimension
        for (int i = 2; i <= nx-1; ++i)
        {
            for (int j = 2; j <= ny-1; ++j)
            {
                for (int k = 2; k <= nz-1; ++k)
                {
                    // calculate the X dim
                    xTerm = Cx * (T[i-1][j][k] + T[i+1][j][k]);
                    // calculate the Y dim
                    yTerm = Cy * (T[i][j-1][k] + T[i][j+1][k]);
                    // calculate the Z dim
                    zTerm = Cz * (T[i][j][k-1] + T[i][j][k+1]);
                    // calculate Tn
                    tnTerm = T[i][j][k] * (1 - (2*Cx + 2*Cy + 2*Cz));
                   
                    if (i == 50 && j == 50 && k == 50)
                    {
                        // combine terms with the Tn term
                        Tnew[i][j][k] = tnTerm + xTerm + yTerm + zTerm;
                        std::cout << "xTerm=" << xTerm << std::endl;
                        std::cout << "yTerm=" << yTerm << std::endl;
                        std::cout << "zTerm=" << zTerm << std::endl;
                        std::cout << "tnTerm=" << tnTerm << std::endl;
                    }
                    //exit(1);
                }
            }
        }
        std::cout << Tnew[50][50][50] << std::endl;

        // set the BCs
        resetBC();
        
        

    }

    void buildLinearSystem(){}
    
private:
    // previous/current timestep
    double***   T;     
    double***   Tnew;
    
    // boundary conditions (sides of a cube)
    double**  xy_0;
    double**  xz_0;
    double**  yz_0;
    double**  xy_N;
    double**  xz_N;
    double**  yz_N;

    // Constants and intervals for calc
    int nx;
    int ny;
    int nz;
    int l_x;
    int l_y;
    int l_z;
    double dx;
    double dy;
    double dz;
    double dt;
    double alpha;
    double Cx;
    double Cy;
    double Cz;
};



