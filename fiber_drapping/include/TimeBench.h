#pragma once
#include <omp.h>
#include <iostream>

class TimeBench
{
	double timeStart;
	double timeEnd;
public:
	TimeBench() { timeStart = omp_get_wtime(); };
	~TimeBench()
	{ 
		timeEnd = omp_get_wtime();
		std::cout << "Time mesurment result: " << timeEnd - timeStart << "\n";
	};
};