/*
* Copyright (C) 2016 by David Hoksza (david.hoksza@gmail.com)
*
* Released under the MIT license, see LICENSE.txt
*/

#ifndef STATS_H
#define STATS_H

#include <cmath>

#include "constants.h"
#include "common.h"

const double PI = 3.141592653589793;
const double E = std::exp(1.0);
const float inv_sqrt_2pi = 0.3989422804014327;
const int max__reasonable_stddev = 20;

#define M

namespace stats {

	SCORE_TYPE log_tab[CNT_PROB_BINS + 1];
	SCORE_TYPE gaussian_map[(2 * max__reasonable_stddev) * CNT_PROB_BINS + 1];
	SCORE_TYPE poisson_map[MAX_OPT_MAP_WINDOW + 1][MAX_FRAGMENT_LENGTH + 1];

	inline int factorial(int n)
	{
		switch (n)
		{
			case 0: return 1;
			case 1: return 1;
			case 2: return 2;
			case 3: return 6;
			case 4: return 24;
			case 5: return 120;
			case 6: return 720;
			case 7: return 5040;
			case 8: return 40320;
			case 9: return 362880;
			case 10: return 3628800;	
			default:
				break;
		}

		//return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
		
		error_exit("Factorial computation only up to 10 available.");
		return 0; //just so that the compiler does not complain
	}

	inline float pdf_gaussian_full(float x, float mean, float stddev)
	{
		//return (1 / (stddev * sqrt(2 * PI))) * exp(-0.5 * pow((x - mean) / stddev, 2.0));
		float a = (x - mean) / stddev;
		return  inv_sqrt_2pi / stddev * exp(-0.5f * a * a);
	}

	float pdf_gaussian(float x, float mean, float stddev)
	{
		if (x < -max__reasonable_stddev || x > max__reasonable_stddev) return 0;
		else{
			float a = gaussian_map[(int)(x * CNT_PROB_BINS + CNT_PROB_BINS * max__reasonable_stddev)];
			return a;
		}
	}

	inline SCORE_TYPE pdf_poisson(int n, int lambda_times_false_cut_p)
	{
		return poisson_map[n][lambda_times_false_cut_p];
	}

	inline float pdf_poisson_full(int n, float lambda)
	{
		return pow(lambda, n)*pow(E, -lambda) / factorial(n);
	}

	SCORE_TYPE transform_prob(SCORE_TYPE p)
	{
		SCORE_TYPE a = log_tab[(int)(p*CNT_PROB_BINS)];
		return a;
	}

	void init_transform_prob()
	{
		for (int p = CNT_PROB_BINS; p > 0; p--) {
			SCORE_TYPE aux = -log((float)p / CNT_PROB_BINS);
			log_tab[p] = (aux > SUB_MAX) ? SUB_MAX : aux;
		}
		log_tab[0] = SUB_MAX;
	}

	void precompute_gaussian()
	{
		for (int i = -CNT_PROB_BINS * max__reasonable_stddev; i <= CNT_PROB_BINS * max__reasonable_stddev; i++)
		{
			gaussian_map[i + CNT_PROB_BINS * max__reasonable_stddev] = pdf_gaussian_full(i / (SCORE_TYPE)CNT_PROB_BINS, 0, 1);
		}

	}

	void precompute_poisson(float lambda_scale)
	{
		for (int i = 0; i <= MAX_OPT_MAP_WINDOW; i++)
		{
			for (int j = 0; j <= MAX_FRAGMENT_LENGTH; j++)
			{
				poisson_map[i][j] = pdf_poisson_full(i, j*lambda_scale);
			}
		}
	}

	void init_stats(float poisson_lambda_scale)
	{
		init_transform_prob();
		precompute_gaussian();
		precompute_poisson(poisson_lambda_scale);
	}
}


#endif // STATS_H