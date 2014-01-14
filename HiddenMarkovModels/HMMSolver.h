#pragma once
class HMMSolver
{
public:
	HMMSolver();
	double get_prob(const ObsSeq &, const Model &) const;
	//void get_model(ObsSeq, int, std::vector<OBS_TYPE>);
	Model get_model(ObsSeq, Model);
	~HMMSolver();
private:
	Model * _m;
	ProbDistrList get_alpha(const ObsSeq &, const Model &, bool) const;
};

