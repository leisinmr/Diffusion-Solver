
#include <math.h>

class Vector
{
public:
    Vector(int size) : m_size(size)
    {
        m_data = (double *)malloc(size*sizeof(double *));
    }
    ~Vector()
    {
        free(m_data);
    }

    double& operator()(int i) 
    {
        return m_data[i];
    }

    double* __internal_state() 
    {
        return m_data;
    }

private:
    int     m_size;
    double* m_data;
};


class Matrix3D
{
public:
    Matrix3D(int x, int y, int z) : m_x(x), m_y(y), m_z(z)
    {
        m_data = (double***)malloc(m_x * sizeof(double **));
        for(int i=0;i<m_x;i++)
        {
            m_data[i] = (double**)malloc(m_y * sizeof(double *));
            for(int j=0;j<m_y;j++)
            {
                m_data[i][j] = (double*)malloc(m_z*sizeof(double));
            }
        }
    }

    ~Matrix3D()
    {
        for(int i=0;i<m_x;i++)
        {
            for(int j=0;j<m_y;j++)
            {
                free(m_data[i][j]);
            }
            free(m_data[i]);
        }
        free(m_data);
    }

    double& operator()(int x, int y, int z)
    {
        return m_data[x][y][z];
    }

    int firstDimSize() const { return m_x; }
    int secondDimSize() const { return m_x; }
    int thirdDimSize() const { return m_x; }

    /*class R
    {
         private:
             friend class M; // Only M can create these objects.
             R(M& parent,int row): m_parent(parent),m_row(row) {}
         public:
              int& operator[](int col) {return m_parent.at(m_row,col);}
         private:
              M&  m_parent;
              int m_row;
    };

    R operator[](int row) {return R(*this,row);}*/
    double*** __internal_state() 
    {
        return m_data;
    }

private:
    int         m_x;
    int         m_y;
    int         m_z;
    double***   m_data;
};

class LinespaceFactory
{
};

void createLinespace(Vector &v, int start, int end, int numIntervals)
{
    double increment = ((double)end - start) / numIntervals;
    for (int i = 0; i < numIntervals; ++i)
        v(i) = increment*i;
}


class GaussianFactory
{
public:
    void operator()(Matrix3D &m, double mu, Vector &x, Vector& y, Vector &z)
    {
        for (int i = 0; i < m.firstDimSize(); ++i)
        {
            for (int j = 0; j < m.secondDimSize(); ++j)
            {
                for (int k = 0; k < m.thirdDimSize(); ++k)
                {
                    m(i,j,k) = exp(-1*pow(5*x(i) - 2.5, 2)) *
                               exp(-1*pow(5*y(i) - 2.5, 2)) *
                               exp(-1*pow(5*z(i) - 2.5, 2));
                }
            }
        }
    }
};

template<typename solution_type, typename BC_type, typename vector_type>
class Heat3D
{
    typedef solution_type Matrix;
    typedef BC_type BC;
    typedef vector_type Vector;

public:

    Heat3D(int _nx, int _ny, int _nz, int _l_x, int _l_y, int _l_z, 
            double _dx, double _dy, double _dz, double _dt, double alpha) 
        : T(_nx, _ny, _nz), Tnew(_nx, _ny, _nz), ic(_nx, _ny, _nz), 
          nx(_nx), ny(_ny), nz(_nz), l_x(_l_x), l_y(_l_y), l_z(_l_z), 
          dx(_dx), dy(_dy), dz(_dz), dt(_dt)
    {
    }  

    void setIC(Matrix ic)
    {
    }

    void setIC(Matrix ic, int xSize, int ySize, int zSize)
    {
    }

    void setBC(BC xy_0, BC xz_0, BC yz_0, BC xy_N, BC xz_N, BC yz_N)
    {
    }

    void setBC(BC xy_0, BC xz_0, BC yz_0, BC xy_N, BC xz_N, BC yz_N, 
               int xSize, int ySize, int zSize)
    {
    }

    void solve(int steps)
    {
    }

private:
    void solveStep()
    {
    }
    
private:
    // previous/current timestep
    //double***   T;     
    //double***   Tnew;
    Matrix    T;
    Matrix    Tnew;
    
    // initial conditions
    Matrix   ic;
   
    // boundary conditions (sides of a cube)
    BC  xy_0;
    BC  xz_0;
    BC  yz_0;
    BC  xy_N;
    BC  xz_N;
    BC  yz_N;

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
    double C;
};



