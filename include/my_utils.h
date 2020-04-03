#ifndef MY_UTILS_H
#define MY_UTILS_H

namespace ACTIONet {
	mat sampleUnif(int l, int m, double a, double b, int seed);
	void gram_schmidt(mat& A);
	field<mat> eigSVD(mat A);
	
	template<class Function>
	void ParallelFor(size_t start, size_t end, size_t numThreads, Function fn);	
}

#endif
