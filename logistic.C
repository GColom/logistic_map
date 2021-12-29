#include<iostream>
#include<vector>
#include"matplotlibcpp.h"

namespace plt = matplotlibcpp;

double * logistic_map_sample(const double & r, const double & x0, const int & niter)
{
	if( x0 > 1. || x0 < 0.)
	{
		std::cout << "Non physical initial coordinate.\n Interrupting!\n" << std::endl;
		return 0;
	}

	double * orbit = new double[niter];
	
	orbit[0] = x0;

	for(int t = 1; t < niter; t++)
	{
		orbit[t] = r * orbit[t - 1 ] * (1 - orbit[t - 1]);
	}

	return orbit;
}

int main(int argc, char const *argv[])
{
	double rr  ;
	double lr  ;
	double dper;
	double init;
	int    N   ;
	int    Nr  ;

	if(argc == 1)
	{
		std::cout << "No input parameters passed, passing default values." << std::endl;
		init = .5;
		N    = 500;
		dper = .925;
		Nr   = 10000;
		lr   = 2.5;
		rr   = 4;
	}

	if(argc == 6)
	{
		rr   = std::stod(argv[6]);
		lr   = std::stod(argv[5]);
		Nr      = std::stoi(argv[4]);
		dper = std::stod(argv[3]);
		N       = std::stoi(argv[2]);
		init = std::stod(argv[1]);
	}
	
	if(argc < 6 && argc > 1)
	{
		std::cout << "Calling signature: "<< argv[0] <<"x0 xniter discard% rniter lr rr"<< std::endl;

		return -1;
	}



	std::cout << "Logistic map orbit diagram for "<< Nr <<" values of r between "<< lr << " and " << rr << ", " << std::endl;
	std::cout << N << " map iterations per value."<<std::endl;
	std::cout << "Discard initial " << dper * 100 << " % of each trajectory to remove transients." << std::endl;

	double ** store = new double * [Nr];

	int t_min = int(dper * N);
	double r_step = (rr-lr)/Nr;
	for(int ir = 0; ir < Nr; ++ir)
	{
		store[ir] = logistic_map_sample(lr + ir * r_step, init, N);
	}

	
	std::vector<double> r(Nr * N,0), x(Nr * N,0);
	for(int i = 0; i < Nr; i++) 
		{
			for(int t = t_min; t < N; t++) 
				{
					r[i * N + t] = lr + i * r_step;
					x[i * N + t] = store[i][t];
				}
		}

	
	plt::figure(1);
	plt::plot(r, x, ".");
	plt::xlabel("r");
	plt::xlim(lr,rr);
	plt::ylabel("x");
	plt::grid(true);
	plt::show();
	/*
	plt::figure(2);
	plt::hist(x, 100);
	plt::grid(true);
	plt::show();
	*/
	return 0;
}