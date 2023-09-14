#ifndef MYTIME_H
#define MYTIME_H

//************************************************************************
// Include
//************************************************************************
#include <signal.h>
#include <iostream.h>
#include <time.h>

//************************************************************************
// Timer
//************************************************************************
class Timer
{
    double dist_time_cpu; 
    double dist_time; 

    clock_t clock_start, clock_end;
    time_t time_start, time_end;
	struct tm *tmstart, *tmend;

public:
    void start();
    void end();
};

//************************************************************************
// Typedef
//************************************************************************
typedef class Timer MYtime;

//************************************************************************
// Function
//************************************************************************
void Timer::start()
{
    clock_start = clock();
    time (&time_start);
	tmstart=localtime(&time_start);
	cout << "start time: "<<asctime(tmstart)<<endl;
}


void Timer::end()
{
    time (&time_end);
    clock_end = clock ();
    dist_time_cpu = (double) (clock_end - clock_start) / (double) CLOCKS_PER_SEC;
	tmend=localtime(&time_end);
    dist_time = difftime(time_end, time_start);
    
	cout <<endl<<endl;
    printf("%f sec cpu,   %f sec real\n", 
	   dist_time_cpu, dist_time);
//	tmstart=localtime(&time_start);
//	cout << "start time: "<<asctime(tmstart);
	cout<< "end time: "<<asctime(tmend)<<endl;
}
#endif