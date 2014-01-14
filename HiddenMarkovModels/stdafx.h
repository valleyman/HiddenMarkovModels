// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#define OBS_TYPE int

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assert.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>

#define TRULY_RDM 0

namespace ublas = boost::numeric::ublas;

typedef std::vector<std::pair<int, OBS_TYPE> > StateObsSeq;
typedef std::vector<OBS_TYPE> ObsSeq;
typedef ublas::vector<double> ProbDistr;
typedef std::vector<ublas::vector<double> > ProbDistrList;

// TODO: reference additional headers your program requires here
