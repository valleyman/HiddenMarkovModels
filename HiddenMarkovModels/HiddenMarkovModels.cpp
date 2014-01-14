// HiddenMarkovModels.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Observable.h"
#include "Model.h"
#include "HMMSolver.h"
#include "Utils.h"
using namespace std;

std::string to_str(matrix<double> &m)
{
	std::string s = "";
	for (int i = 0; i < m.size1(); i++)
	{
		for (int j = 0; j < m.size2(); j++)
		{
			s += m(i, j);
			s += ",";
		}
		s += ";";
	}
	return s;
}

int _tmain(int argc, _TCHAR* argv[])
{

	int endval;

	matrix<double> trans(3, 3, 0.);
	trans(0, 0) = 0.;
	trans(0, 1) = 0.1;
	trans(0, 2) = 1 - trans(0, 0) - trans(0, 1);
	trans(1, 0) = 0.2;
	trans(1, 1) = 0.;
	trans(1, 2) = 1 - trans(1, 0) - trans(1, 1);
	trans(2, 0) = 0.3;
	trans(2, 1) = 0.7;
	trans(2, 2) = 1 - trans(2, 0) - trans(2, 1);

	matrix<double> o_probs(3, 3, 0.);
	o_probs(0, 0) = 0.4;
	o_probs(0, 1) = 0.3;
	o_probs(0, 2) = 1 - o_probs(0, 0) - o_probs(0, 1);
	o_probs(1, 0) = 0.1;
	o_probs(1, 1) = 0.1;
	o_probs(1, 2) = 1 - o_probs(1, 0) - o_probs(1, 1);
	o_probs(2, 0) = 0.1;
	o_probs(2, 1) = 0.1;
	o_probs(2, 2) = 1 - o_probs(2, 0) - o_probs(2, 1);

	ProbDistr init_prob(3);
	init_prob(0) = 0.3;
	init_prob(1) = 0.3;
	init_prob(2) = 1 - init_prob(0) - init_prob(1);

	std::vector<int> abc(3);
	iota(begin(abc), end(abc), 0);


	matrix<double> init_trans(3, 3, 0.);
	init_trans(0, 0) = 0.3;
	init_trans(0, 1) = 0.4;
	init_trans(0, 2) = 1 - init_trans(0, 0) - init_trans(0, 1);
	init_trans(1, 0) = 0.1;
	init_trans(1, 1) = 0.6;
	init_trans(1, 2) = 1 - init_trans(1, 0) - init_trans(1, 1);
	init_trans(2, 0) = 0.1;
	init_trans(2, 1) = 0.1;
	init_trans(2, 2) = 1 - init_trans(2, 0) - init_trans(2, 1);

	matrix<double> init_obs_prob(3, 3, 0.);
	init_obs_prob(0, 0) = 0.2;
	init_obs_prob(0, 1) = 0.6;
	init_obs_prob(0, 2) = 1 - init_obs_prob(0, 0) - init_obs_prob(0, 1);
	init_obs_prob(1, 0) = 0.7;
	init_obs_prob(1, 1) = 0.2;
	init_obs_prob(1, 2) = 1 - init_obs_prob(1, 0) - init_obs_prob(1, 1);
	init_obs_prob(2, 0) = 0.3;
	init_obs_prob(2, 1) = 0.1;
	init_obs_prob(2, 2) = 1 - init_obs_prob(2, 0) - init_obs_prob(2, 1);

	ProbDistr init_init_prob(3);
	init_init_prob(0) = 0.2;
	init_init_prob(1) = 0.4;
	init_init_prob(2) = 1 - init_init_prob(0) - init_init_prob(1);

	HMMSolver s;

	Model m(trans, o_probs, init_prob, abc);
	Model m2(init_trans, init_obs_prob, init_init_prob, abc);

	
	
	for (int iter = 0; iter <= 10; iter++)
	{
		StateObsSeq seq = m.gen_observations_seq(100);

		ObsSeq obss;
		for (std::pair<int, OBS_TYPE> o : seq)
		{
			obss.push_back(o.second);
		}
		cout << "prob m: " << s.get_prob(obss, m) << endl;
		cout << "prob m2: " << s.get_prob(obss, m2) << endl;
		m2 = Model(s.get_model(obss, m2));
		cout << "prob m3: " << s.get_prob(obss, m2) << endl;
		cout << endl << endl;
	}
	
	cout << "Matrice di Transizione" << endl << endl;
	cout << m2.get_trans();
	/*cout << "Sequenza:" << endl;
	for (StateObsSeq::const_iterator it = seq.cbegin(); it != seq.cend(); it++)
	{
		cout << "T: " << std::distance(seq.cbegin(), it) << " | S: " << (*it).first << " | O: " << (*it).second << endl;
	}*/
	
	
	cin >> endval;
	return 0;
}


 
