#ifndef MY_UTILS_H
#define MY_UTILS_H

namespace ACTIONet {
	mat sampleUnif(int l, int m, double a, double b, int seed);
	void gram_schmidt(mat& A);
	field<mat> eigSVD(mat A);
	mat zscore(mat A);
}

#endif
