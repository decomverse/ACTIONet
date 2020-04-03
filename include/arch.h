#ifndef ARCH_H
#define ARCH_H

#include <utils.h>
#include <lsqsplx.h>
#include <projsplx.h>


#define NEW_VERSION

/* **************************
 * Alternating Archetypal Analysis 
 * **************************/

/// Alternating Minimization 
/// Each sub-quadratic programming is solved by ActiveSet Method


void arch_dense(const SPAMS_Matrix<double>& X, const SPAMS_Matrix<double>& Z0, SPAMS_Matrix<double>& Z,  SPAMS_Matrix<double>& A, SPAMS_Matrix<double>& B, const int I1, const int I2, const double lambda2, const double epsilon, const bool computeXtX);

template <typename T>
void arch(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z,  SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const int I1 = 3, const int I2 = 20, const T lambda2 = T(10e-5), const T epsilon = T(10e-5),const bool computeZtZ = true);

template <typename T>
void archRobust(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z,  SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const int I1 = 3, const int I2 = 20, const T lambda2 = T(10e-5), const T epsilon = T(10e-5), const T epsilon2 = T(10e-3),const bool computeZtZ = true);

/// General functions including previous ones. Less parameters and simple use, for Python and Matlab interface

template <typename T>
void archetypalAnalysis(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z0, SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const bool robust =false, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50, const int numThreads=-1);

template <typename T>
void archetypalAnalysis(const SPAMS_Matrix<T>& X, SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& A, SpSPAMS_Matrix<T>& B, const bool robust = false, const T epsilon2 = T(10e-3), const bool computeXtX = false, const int stepsFISTA = 5, const int stepsAS = 50, const bool randominit = true, const int numThreads=-1);

template <typename T>
void decompSimplex(const SPAMS_Matrix<T>& X, const SPAMS_Matrix<T>& Z, SpSPAMS_Matrix<T>& alpha, const bool computerZtZ = false, const int numThreads=-1); 

void decompSimplex_dense(const SPAMS_Matrix<double>& X, const SPAMS_Matrix<double>& Z, SPAMS_Matrix<double>& alpha, const bool computerZtZ = false, const int numThreads=-1); 

#endif
