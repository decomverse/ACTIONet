#include <stdint.h>
#include <utility>
#include <numpy/npy_common.h>
#include <functional>
#include <arma_wrapper.h>
#include <ACTIONet.h>

using aw::npint;
using aw::npdouble;
using aw::intmat;
using aw::dmat;
using aw::dcube; 
using aw::dvec; 

namespace py = pybind11;
using namespace py::literals;

// Computes reduced kernel matrix for a given (single-cell) profile
//
// @param S Input matrix (sparse)
// @param reduced_dim Dimension of the reduced kernel matrix (default=50)
// @param iters Number of SVD iterations (default=5)
// @param seed Random seed (default=0)
// @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION method (1) is implemented (default=1)
// @param SVD_algorithm SVD algorithm to use. Currently supported methods are Halko (1) and Feng (2) (default=1)
// 
// @return A named list with S_r, V, lambda, and exp_var. \itemize{
// \item S_r: reduced kernel matrix of size reduced_dim x #samples.
// \item V: Associated left singular-vectors (useful for reconstructing discriminative scores for features, such as genes).
// \item lambda, exp_var: Summary statistics of the sigular-values.
// }
// 
// @examples
py::dict reduce_kernel(arma::SpMat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, reduction_algorithm, SVD_algorithm);				
	

    py::dict res;
    
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
	return res;
}


// Computes reduced kernel matrix for a given (single-cell) profile
//
// @param S Input matrix (dense)
// @param reduced_dim Dimension of the reduced kernel matrix (default=50)
// @param iters Number of SVD iterations (default=5)
// @param seed Random seed (default=0)
// @param reduction_algorithm Kernel reduction algorithm. Currently only ACTION method (1) is implemented (default=1)
// @param SVD_algorithm SVD algorithm to use. Currently supported methods are Halko (1) and Feng (2) (default=1)
// 
// @return A named list with S_r, V, lambda, and exp_var. \itemize{
// \item S_r: reduced kernel matrix of size reduced_dim x #samples.
// \item V: Associated left singular-vectors (useful for reconstructing discriminative scores for features, such as genes).
// \item lambda, exp_var: Summary statistics of the sigular-values.
// }
// 
// @examples
py::dict reduce_kernel(arma::Mat<npdouble> &S, int reduced_dim = 50, int iters = 5, int seed = 0, int reduction_algorithm = 1, int SVD_algorithm = 1) {
	
	ACTIONet::ReducedKernel reduction;	
	reduction = ACTIONet::reduce_kernel(S, reduced_dim, iters, seed, reduction_algorithm, SVD_algorithm);				
	

    py::dict res;
    
	res["S_r"] = reduction.S_r;		
	res["V"] = reduction.V;
	res["lambda"] = reduction.lambda;
	res["explained_var"] = reduction.exp_var;	
		
	return res;
}




   
   
   

PYBIND11_MODULE(ACTIONet, m) {
    m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

    )pbdoc";
    	
		m.def("reduce_kernel", py::overload_cast<arma::SpMat<npdouble> &, int, int, int, int, int>(&reduce_kernel), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S")=py::none(), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("reduction_algorithm")=1, py::arg("SVD_algorithm")=1);
		m.def("reduce_kernel", py::overload_cast<arma::Mat<npdouble> &, int, int, int, int, int>(&reduce_kernel), "Computes reduced kernel matrix for a given (single-cell) profile", py::arg("S")=py::none(), py::arg("reduced_dim")=50, py::arg("iters")=5, py::arg("seed")=0, py::arg("reduction_algorithm")=1, py::arg("SVD_algorithm")=1);


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif        
}

