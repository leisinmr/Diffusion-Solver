/*
 * Sample 1D Heat diffusion solver
 *
 *
 */

#include <iostream>
#include <cmath>


class FileWriter
{ 
    public:
        FileWriter(std::string name, char delimiter = ',') 
            : m_name(name), m_delim(delimiter)
        {
        }

        ~FileWriter()
        {
            // TODO
        }

        bool open() 
        {
            m_handle = fopen(m_name.c_str(), "w");
        }

        void writeline(std::string line)
        {
            // TODO
        }

        // write any data type converable to a string
        template <typename T>
        void writeline(T items, size_t num)
        {
        }

    private:
        std::string m_name;
        char m_delim;
        FILE *m_handle;
};

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

    double dx = L/nx;
    double dt = 0.0005;

    int nsteps = 1000;

    double C = alpha*dt/pow(dx,2);

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

    for (int n = 0; n < nsteps; ++n)
    {
        for (int i = 2; i < nx+1; ++i)
        {
            Tnew[i] = T[i] + C*(T[i-1] - 2*T[i] + T[i+1]);
        }
    }

}

