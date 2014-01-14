#include "stdafx.h"
#include "Observable.h"
#include "Model.h"

using namespace boost::numeric::ublas;
using namespace boost::assign;

Model::Model(matrix<double> trans, matrix<double> obs_prob, vector<double> init_prob, std::vector<OBS_TYPE> alphabeth) : _obs(obs_prob, alphabeth)//, _trans_gen(trans.size1())
{
	if (trans.size1() != trans.size2() || trans.size1() != obs_prob.size1() || init_prob.size() != trans.size1())
		std::cout << "Inconsistent model parameters!";
		//TODO: Throw Exception
	else
	{
		if (TRULY_RDM == 1) _rng.seed(static_cast<unsigned int>(std::time(0)+1));
		_trans = trans;
		_init_prob = init_prob;
		_states_n = trans.size1();
		_alphabeth = alphabeth;
		_alphabeth_n = alphabeth.size();

		boost::random::discrete_distribution<> start_distr(_init_prob.begin(), _init_prob.end());
		_start_gen = new boost::variate_generator< boost::mt19937, boost::random::discrete_distribution<> >(_rng, start_distr);

		for (unsigned int i = 0; i < _states_n; i++)
		{
			matrix_row<matrix<double> > r((_trans), i);
			boost::random::discrete_distribution<> distr(r.begin(), r.end());
			_trans_gen += boost::variate_generator< boost::mt19937, boost::random::discrete_distribution<> >(_rng, distr);
		}
	}
}

StateObsSeq Model::gen_observations_seq(unsigned int tot_time)
{
	StateObsSeq seq(tot_time);
	
	int state = (*_start_gen)();
	seq[0] = std::pair<int,OBS_TYPE>(state, _obs.measure(state));
	for (unsigned int i = 1; i < tot_time; ++i)
	{
		state = _trans_gen[seq[i-1].first]();
		seq[i] = std::pair<int, OBS_TYPE>(state, _obs.measure(state));
	}
	return seq;
}

// Getters & Setters

Model& Model::set_trans(matrix<double> trans)
{
	_trans = trans;
	return *this;
}
Model& Model::set_trans(int i, int j, double t)
{
	_trans(i,j) = t;
	return *this;
}
matrix<double> Model::get_trans() const
{
	return _trans;
}
double Model::get_trans(int i, int j) const
{
	return _trans(i, j);
}
Model& Model::set_obs_prob(matrix<double> obs_prob)
{
	_obs.set_prob(obs_prob);
	return *this;
}
Model& Model::set_obs_prob(int i, int j, double p)
{
	_obs.set_prob(i, j, p);
	return *this;
}
matrix<double> Model::get_obs_prob() const
{
	return _obs.get_prob();
}
double Model::get_obs_prob(int i, int j) const
{
	return _obs.get_prob(i, j);
}
Model& Model::set_init_prob(vector<double> init_prob)
{
	_init_prob = init_prob;
	return *this;
}
Model& Model::set_init_prob(int i, double p)
{
	_init_prob[i] = p;
	return *this;
}
vector<double> Model::get_init_prob() const
{
	return _init_prob;
}
double Model::get_init_prob(int i) const
{
	return _init_prob[i];
}
int Model::get_states_n() const
{
	return _states_n;
}
std::vector<OBS_TYPE> Model::get_alphabeth() const
{
	return _alphabeth;
}

Model::~Model()
{
}
