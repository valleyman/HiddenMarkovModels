#pragma once
#include <stdio.h>
#include <tchar.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/assign/std/vector.hpp>

namespace ublas = boost::numeric::ublas;

#define MIN_CHOP_VALUE 0.0000001

namespace Utils
{
	template <typename T>
	void dump_std_std(std::vector<std::vector<T> > m)
	{
		std::cout << "Matrice" << std::endl;
		std::cout << std::endl;
		for (std::vector<T> r : m)
		{
			for (T e : r)
			{
				std::cout << e << "\t";
			}
			std::cout << "\n";
		}
	}

	template <typename T>
	void dump_std_bst(std::vector<boost::numeric::ublas::vector<T> > m)
	{
		std::cout << "Matrice" << std::endl;
		std::cout << std::endl;
		for (ublas::vector<T> r : m)
		{
			for (T e : r)
			{
				std::cout << e << "\t";
			}
			std::cout << "\n";
		}
	}

	void chop(boost::numeric::ublas::matrix<double> &m)
	{
		for (double& e : m.data())
		{
			if (abs(e) <= MIN_CHOP_VALUE) e = 0.;
		}
	}
	void chop(std::vector<double> &v)
	{
		for (double& e : v)
		{
			if(abs(e) <= MIN_CHOP_VALUE) e = 0.;
		}
	}
	double sum(ublas::vector<double> v)
	{
		double sum = 0;
		for (int i = 0; i < v.size(); i++) sum += v[i];
		return sum;
	}
}