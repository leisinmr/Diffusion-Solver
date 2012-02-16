#include <sys/time.h>

class Timer
{
    typedef struct timeval time_type;
public:
    Timer() : tv1(), tv2(), microseconds(0.0)
    {
    }

    ~Timer() {}

    void start()
    {
        gettimeofday(&tv1, NULL);
    }

    double stop()
    {
        gettimeofday(&tv2, NULL);
        microseconds = (double)(tv2.tv_sec - tv1.tv_sec) * 1000; // sec to ms
        microseconds += (double)(tv2.tv_usec - tv1.tv_usec) / 1000; // ms to us

        return microseconds;
    }

    double elapsedTime()
    {
        return microseconds;
    }


private:
    time_type tv1;
    time_type tv2;
    double microseconds;
};


