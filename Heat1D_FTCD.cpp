/*
 * Sample 1D Heat diffusion solver
 * 
 * For testing initial and boundary conditions against Matlab
 *
 */

#include <iostream>
#include <cmath>
#include <fstream>

class FileStream
{ 
    public:
        FileStream(std::string name, char delimiter = ',') 
            : m_name(name), m_delim(delimiter), m_stream(name.c_str())
        {
        }

        ~FileStream()
        {
            m_stream.close();
        }

        void writeline(double *row, int size)
        {
            for (int i = 0; i < size; ++i)
            {
                m_stream << row[i] << ",";
                //m_stream << (i != size-1) ? "," : "";
            }
            m_stream << '\n';
        }

        // write any data type converable to a string
        template <typename T>
        void writeline(T items, size_t num)
        {

        }

    private:
        std::string m_name;
        char m_delim;
        std::ofstream m_stream;
};

/*template <typename stream>
stream& operator<<(const stream& s, double val)
{
    stream << val;
}*/

double* createLinspace(int start, int end, int numIntervals)
{
    double increment = ((double)end - start) / numIntervals;
    double *linspace = (double *)malloc(sizeof(double)*numIntervals);

    if (linspace == NULL)
        return NULL;

    for (int i = 0; i < numIntervals; ++i) 
    {
        linspace[i] = increment*i;
    }

    return linspace;
}


double* createGaussian(double mu, double *x, int sizeX)
{
    double *gaussian = (double *)malloc(sizeof(double)*sizeX);

    if (gaussian == NULL)
        return NULL;

    for (int i = 0; i < sizeX; ++i)
    {
        gaussian[i] = exp(-1*pow(5*x[i] - 2.5, 2));
    }

    return gaussian;
}

int main(int argc, char** argv)
{
    int nx = 1000;
    int L = 1;
    double alpha = 0.001;

    double dx = (double)L/nx;
    double dt = 0.0005;

    int nsteps = 1000;

    double C = alpha*dt/pow(dx,2);

    std::cout << "Starting 1D Heat Solver with: " << std::endl;
    std::cout << "nx = " << nx << ", alpha= " << alpha;
    std::cout << ", dx= " << dx << ", dt= " << dt << std::endl;

    std::cout << "C = " << C << std::endl;

    double *x = createLinspace(0, L, nx+2);

    if (x == NULL) 
    {
        std::cout << "Problem creating linspace" << std::endl;
        exit(2);
    }
  
    double *T = createGaussian(0.0, x, nx+2);
    double *Tnew = (double *)malloc(sizeof(double)*(nx+2));

    if (T == NULL)
    {
        std::cout << "Problem creating gaussian intial condition" << std::endl;
        exit(2);
    }

    FileStream f("Heat1D.csv");
    for (int n = 0; n < nsteps; ++n)
    {
        // Zero at the boundary conditions
        Tnew[0] = 0; Tnew[nx+1] = 0;
        // solve for one temperature step
        for (int i = 2; i < nx+1; ++i)
        {
            // apply the FTCD across the 1D rod
            Tnew[i] = T[i] + C*(T[i-1] - 2*T[i] + T[i+1]);

        }

        // write out to a file
        f.writeline(Tnew, nx+2);
        
        for (int i = 2; i < nx+1; ++i)
        {
            // set Tnew to T so we can do the next step
            T[i] = Tnew[i];
        }
    }

}

