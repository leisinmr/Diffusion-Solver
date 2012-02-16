#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <fstream>

#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*

template<typename T>
T*** allocate3DMatrix(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	T ***t;

	/* allocate pointers to pointers to rows */
	t=(T ***) malloc((size_t)((nrow+NR_END)*sizeof(T**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(T **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(T*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(T *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(T)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

template<typename T>
void free_3tensor(T ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

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

        // write csv output
        void writeline(double *row, int size)
        {
            for (int i = 0; i < size; ++i)
            {
                m_stream << row[i] << ",";
            }
            m_stream << '\n';
        }

        // direct access
        std::ofstream& _stream() { return m_stream; }

    private:
        std::string m_name;
        char m_delim;
        std::ofstream m_stream;
};


#endif // __UTIL_U__

