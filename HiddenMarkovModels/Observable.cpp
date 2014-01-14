#include "stdafx.h"
#include "Observable.h"

using namespace boost::assign;

Observable::Observable(const matrix<double> &prob, const std::vector<OBS_TYPE> &alphabeth) //: _obs_gen(prob.size1())
{
	_prob = prob;
	_alphabeth = alphabeth;
	if (TRULY_RDM == 1) _rng.seed(static_cast<unsigned int>(std::time(0)));
	for (unsigned int i = 0; i < prob.size1(); i++)
	{
		matrix_row<matrix<double> > r((_prob), i);
		boost::random::discrete_distribution<> distr(r.begin(), r.end());
		_obs_gen += boost::variate_generator< boost::mt19937, boost::random::discrete_distribution<> >(_rng, distr);
	}
}

OBS_TYPE Observable::measure(unsigned int state)
{
	if (state < _obs_gen.size())
		return _alphabeth[_obs_gen[state]()];
	else
		// TODO: Throw exception
		return NULL;
}

// Getters & Setters

Observable& Observable::set_prob(matrix<double> prob)
{
	_prob = prob;
	return *this;
}
Observable& Observable::set_prob(int i, int j, double p)
{
	_prob(i, j) = p;
	return *this;
}
matrix<double> Observable::get_prob() const
{
	return _prob;
}
double Observable::get_prob(int i, int j) const
{
	return _prob(i, j);
}
int Observable::get_symb_ind(OBS_TYPE symbol) const
{
	auto ind = std::find(_alphabeth.begin(), _alphabeth.end(), symbol);
	return std::distance(_alphabeth.begin(), ind);
}

double Observable::get_symb_prob(int state, OBS_TYPE symbol) const
{
	return _prob(state, get_symb_ind(symbol));
}

Observable::~Observable()
{
}
















