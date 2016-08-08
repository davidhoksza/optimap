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


namespace stats {

	const double PI = 3.141592653589793;
	const double E = std::exp(1.0);
	const float inv_sqrt_2pi = 0.3989422804014327;
	const float sqrt_2pi = 2.50662827463;
	const int max__reasonable_stddev = 20;
	const int LAPLACE_GRANULARITY = 10000;
	const int LAPLACE_MAX = 3; //used in sizing error (ref/exp) -> 3x difference in exp vs ref is considered impossible (stddev is about 0.1)
		


	SCORE_TYPE log_tab[CNT_PROB_BINS + 1];
	SCORE_TYPE gaussian_map[(2 * max__reasonable_stddev) * CNT_PROB_BINS + 1];
	SCORE_TYPE poisson_map[MAX_OPT_MAP_WINDOW + 1][MAX_FRAGMENT_LENGTH + 1];
	SCORE_TYPE poisson_map_0[MAX_OPT_MAP_WINDOW + 1][MAX_FRAGMENT_LENGTH + 1];
	SCORE_TYPE poisson_map_1[MAX_OPT_MAP_WINDOW + 1][MAX_FRAGMENT_LENGTH + 1];
	SCORE_TYPE poisson_map_3[MAX_OPT_MAP_WINDOW + 1][MAX_FRAGMENT_LENGTH + 1];
	SCORE_TYPE laplace_map[4][2 * LAPLACE_MAX * LAPLACE_GRANULARITY + 1]; //laplace distribution for fragments < 2.4 kb (location = 0.858181;	scale = 0.180196)
	//SCORE_TYPE laplace_map_lt_3600[2 * LAPLACE_MAX * LAPLACE_GRANULARITY + 1]; //laplace distribution for fragments < 3.6 kb (location = 0.980760;	scale = 0.071176)
	//SCORE_TYPE laplace_map_lt_4800[2 * LAPLACE_MAX * LAPLACE_GRANULARITY + 1]; //laplace distribution for fragments < 4.8 kb (location = 1.003354;	scale = 0.052800)
	//SCORE_TYPE laplace_map_ge_4800[2 * LAPLACE_MAX * LAPLACE_GRANULARITY + 1]; //laplace distribution for fragments >= 4.8 kb (location = 1.00482;	scale = 0.042428)

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

	inline float pdf_gaussian(float x, float mean, float stddev)
	{
		if (x < -max__reasonable_stddev || x > max__reasonable_stddev) return 0;
		else return gaussian_map[(int)(x * CNT_PROB_BINS) + CNT_PROB_BINS * max__reasonable_stddev]; //+ CNT_PROB_BINS * max__reasonable_stddev -> index starts from 0
	}

	inline float pdf_laplace_full(float x, float location, float scale)
	{
		//****************************************************************************80
		//
		//  Purpose:
		//
		//    LAPLACE_PDF evaluates the Laplace PDF.
		//
		//  Discussion:
		//
		//    PDF(A,B;X) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
		//
		//  Discussion:
		//
		//    The Laplace PDF is also known as the Double Exponential PDF.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    09 February 1999
		//
		//  Author:
		//
		//    John Burkardt
		//
		//  Parameters:
		//
		//    Input, double X, the argument of the PDF.
		//
		//    Input, double A, B, the parameters of the PDF.
		//    0.0 < B.
		//
		//    Output, double LAPLACE_PDF, the value of the PDF.
		//
		double pdf;

		pdf = exp(-fabs(x - location) / scale) / (2.0 * scale);

		return pdf;

	}

	inline SCORE_TYPE pdf_poisson(int n, int lambda_times_false_cut_p)
	{	
		if (lambda_times_false_cut_p > MAX_FRAGMENT_LENGTH) {
			std::ostringstream ss;
			ss << "Maximum fragment length currently set to " << MAX_FRAGMENT_LENGTH;
			error_exit(ss.str());
		}
		return poisson_map[n][lambda_times_false_cut_p];
	}

	inline SCORE_TYPE pdf_poisson_200kb(int lambda_200kb, int k, int frag_length)
	{
		if (frag_length > MAX_FRAGMENT_LENGTH) {
			std::ostringstream ss;
			ss << "Maximum fragment length currently set to " << MAX_FRAGMENT_LENGTH;
			error_exit(ss.str());
		}
		switch (lambda_200kb)
		{
		case 0:
			return poisson_map_0[k][frag_length];
		case 1:
			return poisson_map_1[k][frag_length];
			
		case 3:
			return poisson_map_3[k][frag_length];
		default:
			std::ostringstream ss;
			ss << "Unsupported value of lambda (" << lambda_200kb << ") for the Poisson distribution.";
			error_exit(ss.str());
			break;
		}		
	}

	inline SCORE_TYPE pdf_laplace(float x, int type)
	{
		if (x < -LAPLACE_MAX || x > LAPLACE_MAX) return 0;
		
		return laplace_map[type][LAPLACE_MAX * LAPLACE_GRANULARITY + (int)(x * LAPLACE_GRANULARITY)];
	}

	inline float pdf_poisson_full(int n, float lambda)
	{
		return pow(lambda, n)*pow(E, -lambda) / factorial(n);
	}

	SCORE_TYPE transform_prob(SCORE_TYPE p)
	{
		int ix = (int)(p*CNT_PROB_BINS);
		//if (ix > CNT_PROB_BINS + 1)
		SCORE_TYPE a = log_tab[ix];
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

	void precompute_poisson_200kb(int lambda_200kb)
	{	
		for (int i = 0; i <= MAX_OPT_MAP_WINDOW; i++)
		{
			for (int j = 0; j <= MAX_FRAGMENT_LENGTH; j++)
			{
				float p = pdf_poisson_full(i, j * (lambda_200kb / 200000.0));
				switch (lambda_200kb)
				{
				case 0:
					poisson_map_0[i][j] = p;
					break;
				case 1:
					poisson_map_1[i][j] = p;
					break;
				case 3:
					poisson_map_3[i][j] = p;
					break;
				default:
					std::ostringstream ss;
					ss << "Unsupported value of lambda (" << lambda_200kb << ") for the Poisson distribution.";
					error_exit(ss.str());
					break;
				}
			}
		}
	}

	void precompute_laplace()
	{
		float locations[] = { 0.858181, 0.980760, 1.003354, 1.00482 };
		float scales[] = { 0.180196, 0.071176, 0.052800, 0.042428 };
		float step = 1.0 / LAPLACE_GRANULARITY;

		int ix_center = LAPLACE_MAX * LAPLACE_GRANULARITY;
		for (int i = 0; i < sizeof(locations) / sizeof(*locations); i++)
		{
			for (int j = 0; j <= LAPLACE_MAX * LAPLACE_GRANULARITY; j++)
			{
				laplace_map[i][ix_center + j] = pdf_laplace_full(j * step, locations[i], scales[i]);
				laplace_map[i][ix_center - j] = pdf_laplace_full(-j * step, locations[i], scales[i]);
			}
		}
	}

	void init_stats(float poisson_lambda_scale)
	{
		init_transform_prob();
		precompute_gaussian();
		precompute_poisson(poisson_lambda_scale);
		precompute_poisson_200kb(0);
		precompute_poisson_200kb(1);
		precompute_poisson_200kb(3);
		precompute_laplace();
	}
}


#endif // STATS_H