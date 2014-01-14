#pragma once

using namespace ublas;

class Model
{
public:
	Model(matrix<double>, matrix<double>, ProbDistr, std::vector<OBS_TYPE>);

	StateObsSeq gen_observations_seq(unsigned int);

	// Getters & Setters
	Model& set_trans(matrix<double>);
	Model& set_trans(int, int, double);
	matrix<double> get_trans() const;
	double get_trans(int, int) const;
	Model& set_obs_prob(matrix<double>);
	Model& set_obs_prob(int, int, double);
	matrix<double> get_obs_prob() const;
	double get_obs_prob(int, int) const;
	Model& set_init_prob(ProbDistr);
	Model& set_init_prob(int, double);
	ProbDistr get_init_prob() const;
	double get_init_prob(int) const;
	int get_states_n() const;
	std::vector<OBS_TYPE> get_alphabeth() const;

	~Model();

private:
	int _states_n, _alphabeth_n;
	boost::mt19937 _rng;
	std::vector<OBS_TYPE> _alphabeth;
	matrix<double> _trans;
	ProbDistr _init_prob;
	Observable _obs;
	boost::variate_generator<boost::mt19937, boost::random::discrete_distribution<> > * _start_gen; // Pointer just to avoid being forced to 
																									// initialize _start_gen in initialization list
	std::vector<boost::variate_generator<boost::mt19937, boost::random::discrete_distribution<> > > _trans_gen;

};

