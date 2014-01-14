#include "stdafx.h"
#include "Observable.h"
#include "Model.h"
#include "HMMSolver.h"

namespace Utils
{
	double sum(ublas::vector<double> v);
}

HMMSolver::HMMSolver()
{
}

ProbDistrList HMMSolver::get_alpha(const ObsSeq &obs_seq, const Model &model, bool scaled) const
{
	int states_num = model.get_states_n();
	int tot_time = obs_seq.size();
	
	ProbDistr init_prob = model.get_init_prob();
	matrix<double> trans = model.get_trans();
	Observable obs = Observable(model.get_obs_prob(), model.get_alphabeth());

	ProbDistrList alpha(tot_time, ProbDistr(states_num));

	for (int i = 0; i < states_num; ++i)
	{
		alpha[0][i] = init_prob[i] * obs.get_symb_prob(i, obs_seq[0]);
	}
	for (int t = 1; t < tot_time; ++t) 
	{
		for (int j = 0; j < states_num; ++j)
		{
			double sum_alpha = 0;
			for (int i = 0; i < states_num; ++i)
			{
				sum_alpha += alpha[t - 1][i] * trans(i, j);
			}
			alpha[t][j] = obs.get_symb_prob(j, obs_seq[t]) * sum_alpha;

		}
		
	}
	return alpha;
}

double HMMSolver::get_prob(const ObsSeq &obs_seq, const Model &model) const
{
	ProbDistrList alpha = get_alpha(obs_seq, model, true);
	double prob = 0;
	for (int i = 0; i < model.get_states_n(); i++)
	{
		prob += alpha[obs_seq.size()-1][i];
	}
	return prob;
}
/*
void HMMSolver::get_model(ObsSeq obs_seq, int states_num, std::vector<OBS_TYPE> alphabeth)
{
	int tot_time = obs_seq.size();
	int alphabeth_size = alphabeth.size();

	//matrix<double> trans(states_num, states_num, 1. / ((double)states_num));
	//matrix<double> obs_prob(states_num, alphabeth_size, 1. / ((double)alphabeth_size));
	//vector<double> init_prob(states_num, 1. / ((double)states_num));


	matrix<double> trans(3, 3, 0.);
	trans(0, 0) = 0.3;
	trans(0, 1) = 0.4;
	trans(0, 2) = 1 - trans(0, 0) - trans(0, 1);
	trans(1, 0) = 0.1;
	trans(1, 1) = 0.6;
	trans(1, 2) = 1 - trans(1, 0) - trans(1, 1);
	trans(2, 0) = 0.1;
	trans(2, 1) = 0.1;
	trans(2, 2) = 1 - trans(2, 0) - trans(2, 1);

	matrix<double> obs_prob(3, 3, 0.);
	obs_prob(0, 0) = 0.2;
	obs_prob(0, 1) = 0.6;
	obs_prob(0, 2) = 1 - obs_prob(0, 0) - obs_prob(0, 1);
	obs_prob(1, 0) = 0.7;
	obs_prob(1, 1) = 0.2;
	obs_prob(1, 2) = 1 - obs_prob(1, 0) - obs_prob(1, 1);
	obs_prob(2, 0) = 0.3;
	obs_prob(2, 1) = 0.1;
	obs_prob(2, 2) = 1 - obs_prob(2, 0) - obs_prob(2, 1);

	ProbDistr init_prob(3);
	init_prob(0) = 0.2;
	init_prob(1) = 0.4;
	init_prob(2) = 1 - init_prob(0) - init_prob(1);

	
	Observable obs(obs_prob, alphabeth);
	
	Model model(trans, obs_prob, init_prob, alphabeth);
	
	for (int iters = 0; iters <= 10; iters++)
	{
		matrix<double>trans(model.get_trans());
		matrix<double>obs_prob(model.get_obs_prob());
		vector<double>init_prob(model.get_init_prob());
		Observable obs(obs_prob, alphabeth);

		ProbDistrList alpha(tot_time, ProbDistr(states_num));
		ProbDistrList beta(tot_time, ProbDistr(states_num));
		ProbDistrList gamma(tot_time, ProbDistr(states_num));
		std::vector<matrix<double> > xi(tot_time, matrix<double>(states_num, states_num));
		std::vector<double> scaling(tot_time);

		// Alpha and Beta
		for (int i = 0; i < states_num; ++i)
		{
			alpha[0][i] = init_prob[i] * obs.get_symb_prob(i, obs_seq[0]);
		}
		for (int t = 1; t < tot_time; ++t)
		{
			double c = 0;
			for (int j = 0; j < states_num; ++j)
			{
				double sum_alpha = 0;
				for (int i = 0; i < states_num; ++i)
				{
					sum_alpha += alpha[t - 1][i] * trans(i, j);
				}
				alpha[t][j] = obs.get_symb_prob(j, obs_seq[t]) * sum_alpha;
				c += alpha[t][j];
			}
			scaling[t] = c;
			for (int j = 0; j < states_num; ++j)
			{
				//std::cout << "alpha:" << alpha[t][j] << std::endl;
				alpha[t][j] /= c;
				//std::cout << "alphaC:" << alpha[t][j] << std::endl;
			}

		}
		for (int i = 0; i < states_num; ++i)
		{
			beta[tot_time - 1][i] = 1;
		}

		for (int t = tot_time - 2; t >= 0; --t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				double sum_beta = 0;
				for (int j = 0; j < states_num; ++j)
				{
					sum_beta += beta[t + 1][j] * trans(i, j) * obs.get_symb_prob(j, obs_seq[t + 1]);
				}
				beta[t][i] = sum_beta / scaling[t];
			}
		}

		// Gamma and Xi
		for (int t = 0; t < tot_time; ++t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				double sum_gamma = 0;
				for (int j = 0; j < states_num; ++j)
				{
					sum_gamma += alpha[t][j] * beta[t][j];
				}
				gamma[t][i] = alpha[t][i] * beta[t][i] / sum_gamma;
			}
		}

		for (int t = 0; t < tot_time - 1; ++t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				for (int j = 0; j < states_num; ++j)
				{
					double sum_xi_k = 0;
					for (int k = 0; k < states_num; ++k)
					{
						double sum_xi_l = 0;
						for (int l = 0; l < states_num; ++l)
						{
							sum_xi_l += alpha[t][k] * trans(k, l) * beta[t + 1][l] * obs.get_symb_prob(l, obs_seq[t + 1]);
						}
						sum_xi_k += sum_xi_l;
					}
					xi[t](i, j) = alpha[t][i] * trans(i, j) * beta[t + 1][j] * obs.get_symb_prob(j, obs_seq[t + 1]) / sum_xi_k;
				}
			}
		}

		// Model parameters update

		// TODO: The assignment of new variables can be done directly on the model instance, for performance purposes
		//std::cout << gamma[0];
		ProbDistr new_init_prob(gamma[0]);
		matrix<double> new_trans(states_num, states_num);

		for (int i = 0; i < states_num; ++i)
		{
			for (int j = 0; j < states_num; ++j)
			{
				double sum_xi = 0;
				double sum_gamma = 0;
				for (int t = 0; t < tot_time - 1; ++t)
				{
					sum_xi += xi[t](i, j);
					sum_gamma += gamma[t][i];
				}
				new_trans(i, j) = sum_xi / sum_gamma;
			}
		}

		//std::cout << std::endl << std::endl << new_trans <<std::endl <<std::endl;

		matrix<double> new_obs_prob(states_num, alphabeth_size);
		for (int i = 0; i < states_num; ++i)
		{
			for (int k = 0; k < alphabeth_size; ++k)
			{
				// TODO: sum_gamma1[i] = previous sum_gamma[i] loop + gamma[t][i] , consider saving variable value for performance purposes.
				double sum_gamma1 = 0;
				double sum_gamma2 = 0;
				for (int t = 0; t < tot_time; ++t)
				{
					if (obs_seq[t] == alphabeth[k])
					{
						sum_gamma1 += gamma[t][i];
					}
					sum_gamma2 += gamma[t][i];
				}
				new_obs_prob(i, k) = sum_gamma1 / sum_gamma2;
			}
		}
		
		//std::cout << new_trans;
		//std::cout << new_init_prob;
		//std::cout << new_obs_prob;

		model.set_trans(new_trans);
		model.set_init_prob(new_init_prob);
		model.set_obs_prob(new_obs_prob);
	}
	
	std::cout << model.get_init_prob();
	std::cout << model.get_trans();
	std::cout << model.get_obs_prob();

}
*/
Model HMMSolver::get_model(ObsSeq obs_seq, Model initm)
{
	int tot_time = obs_seq.size();
	int alphabeth_size = initm.get_alphabeth().size();
	int states_num = initm.get_states_n();
	/*matrix<double> trans(states_num, states_num, 1. / ((double)states_num));
	matrix<double> obs_prob(states_num, alphabeth_size, 1. / ((double)alphabeth_size));
	vector<double> init_prob(states_num, 1. / ((double)states_num));*/


	/*matrix<double> trans = initm.get_trans();
	ProbDistr init_prob = initm.get_init_prob();
	matrix<double> obs_prob = initm.get_obs_prob();
	std::vector<OBS_TYPE> alphabeth = initm.get_alphabeth();*/

//	Observable obs(obs_prob, alphabeth);

	Model model(initm.get_trans(), initm.get_obs_prob(), initm.get_init_prob(), initm.get_alphabeth());

	for (int iters = 0; iters <= 20; iters++)
	{
		matrix<double>trans(model.get_trans());
		matrix<double>obs_prob(model.get_obs_prob());
		vector<double>init_prob(model.get_init_prob());
		Observable obs(obs_prob, model.get_alphabeth());

		ProbDistrList alpha(tot_time, ProbDistr(states_num));
		ProbDistrList beta(tot_time, ProbDistr(states_num));
		ProbDistrList gamma(tot_time, ProbDistr(states_num));
		std::vector<matrix<double> > xi(tot_time, matrix<double>(states_num, states_num));
		std::vector<double> scaling(tot_time);
		
		
		std::cout << "T: " << trans <<std::endl;
		std::cout << "O: " << obs_prob << std::endl;
		std::cout << "I: " << init_prob << std::endl << std::endl;

		// Alpha and Beta
		double s = 0;
		for (int i = 0; i < states_num; ++i)
		{
			alpha[0][i] = init_prob[i] * obs.get_symb_prob(i, obs_seq[0]);
			s += alpha[0][i];
		}
		scaling[0] = s;
		for (int j = 0; j < states_num; ++j)
		{
			alpha[0][j] /= s;
		}
		for (int t = 1; t < tot_time; ++t)
		{
			s = 0;
			for (int j = 0; j < states_num; ++j)
			{
				double sum_alpha = 0;
				for (int i = 0; i < states_num; ++i)
				{
					sum_alpha += alpha[t - 1][i] * trans(i, j);
				}
				alpha[t][j] = obs.get_symb_prob(j, obs_seq[t]) * sum_alpha;
				s += alpha[t][j];
			}
			scaling[t] = s;
			for (int j = 0; j < states_num; ++j)
			{
				//std::cout << "alpha:" << alpha[t][j] << std::endl;
				alpha[t][j] /= s;
				//std::cout << "alphaC:" << alpha[t][j] << std::endl;
			}

		}
		for (int i = 0; i < states_num; ++i)
		{
			beta[tot_time - 1][i] = 1;
		}

		for (int t = tot_time - 2; t >= 0; --t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				double sum_beta = 0;
				for (int j = 0; j < states_num; ++j)
				{
					sum_beta += beta[t + 1][j] * trans(i, j) * obs.get_symb_prob(j, obs_seq[t + 1]);
				}
				
				beta[t][i] = sum_beta / scaling[t];
				//std::cout << "beta:" << beta[t][i] * scaling[t] << std::endl;
				//std::cout << "betaC:" << beta[t][i] << std::endl;
			}
		}

		// Gamma and Xi
		for (int t = 0; t < tot_time; ++t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				double sum_gamma = 0;
				for (int j = 0; j < states_num; ++j)
				{
					sum_gamma += alpha[t][j] * beta[t][j];
				}
				gamma[t][i] = alpha[t][i] * beta[t][i] / sum_gamma;
			}
		}

		for (int t = 0; t < tot_time - 1; ++t)
		{
			for (int i = 0; i < states_num; ++i)
			{
				for (int j = 0; j < states_num; ++j)
				{
					double sum_xi_k = 0;
					for (int k = 0; k < states_num; ++k)
					{
						double sum_xi_l = 0;
						for (int l = 0; l < states_num; ++l)
						{
							sum_xi_l += alpha[t][k] * trans(k, l) * beta[t + 1][l] * obs.get_symb_prob(l, obs_seq[t + 1]);
						}
						sum_xi_k += sum_xi_l;
					}
					xi[t](i, j) = alpha[t][i] * trans(i, j) * beta[t + 1][j] * obs.get_symb_prob(j, obs_seq[t + 1]) / sum_xi_k;
				}
			}
		}

		// Model parameters update

		// TODO: The assignment of new variables can be done directly on the model instance, for performance purposes
		//std::cout << gamma[0];
		ProbDistr new_init_prob(gamma[0]);
		matrix<double> new_trans(states_num, states_num);

		for (int i = 0; i < states_num; ++i)
		{
			for (int j = 0; j < states_num; ++j)
			{
				double sum_xi = 0;
				double sum_gamma = 0;
				for (int t = 0; t < tot_time - 1; ++t)
				{
					sum_xi += xi[t](i, j);
					sum_gamma += gamma[t][i];
				}
				new_trans(i, j) = sum_xi / sum_gamma;
			}
		}

		//std::cout << std::endl << std::endl << new_trans << std::endl << std::endl;

		matrix<double> new_obs_prob(states_num, alphabeth_size);
		for (int i = 0; i < states_num; ++i)
		{
			for (int k = 0; k < alphabeth_size; ++k)
			{
				// TODO: sum_gamma1[i] = previous sum_gamma[i] loop + gamma[t][i] , consider saving variable value for performance purposes.
				double sum_gamma1 = 0;
				double sum_gamma2 = 0;
				for (int t = 0; t < tot_time; ++t)
				{
					if (obs_seq[t] == model.get_alphabeth()[k])
					{
						sum_gamma1 += gamma[t][i];
					}
					sum_gamma2 += gamma[t][i];
				}
				new_obs_prob(i, k) = sum_gamma1 / sum_gamma2;
			}
		}

		//std::cout << new_trans;
		//std::cout << new_init_prob;
		//std::cout << new_obs_prob;

		/*Utils::chop(new_trans);
		Utils::chop(new_init_prob);
		Utils::chop(new_obs_prob);*/

		model.set_trans(new_trans);
		model.set_init_prob(new_init_prob);
		model.set_obs_prob(new_obs_prob);
	}

	return model;

}

HMMSolver::~HMMSolver()
{
}
