#pragma once

using namespace boost::numeric::ublas;
using namespace boost::assign;

class Observable
{
public:
	Observable(const matrix<double> &, const std::vector<OBS_TYPE> &);
	OBS_TYPE measure(unsigned int);
	Observable& set_prob(matrix<double>);
	Observable& set_prob(int, int, double);
	matrix<double> get_prob() const;
	double get_prob(int, int) const;
	int get_symb_ind(OBS_TYPE) const;
	double get_symb_prob(int, OBS_TYPE) const;
	~Observable();
private:
	matrix<double> _prob;
	std::vector<OBS_TYPE> _alphabeth;
	boost::mt19937 _rng;
	std::vector< boost::variate_generator< boost::mt19937, boost::random::discrete_distribution<> > > _obs_gen;
};

