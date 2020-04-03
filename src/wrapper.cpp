#include <stdint.h>
#include <utility>
#include <numpy/npy_common.h>
#include <functional>
#include <arma_wrapper.h>

using aw::npint;
using aw::npdouble;
using aw::intmat;
using aw::dmat;
using aw::dcube; 
using aw::dvec; 



PYBIND11_MODULE(ACTIONet, m) {
    m.doc() = R"pbdoc(
        ACTIONet package
        -----------------------

        .. currentmodule:: ACTIONet

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";
    	
        m.def("sum", [] (const dmat & arr) {	return arma::accu(arr);},
               "Sums the elements in the array.");

        m.def("expmat", [](const dmat& arg) {return dmat(arma::expmat(arg));}); 

        m.def("identity", [](const arma::SpMat<npdouble>& arg) { return arg;});

        m.def("col_identity", [] (const dvec& arg) {return arg;});
    
        m.def("row_identity", [](const arma::Row<npdouble>& arg) {return arg;});

        m.def("cube_identity", [] (const dcube& arg) {return arg;});

        m.def("mat_identity", [] (const dmat& arg) {return arg;});
        m.def("intmat_identity", [] (const intmat& arg) {return arg;});

        m.def("duplicate", [](const arma::SpMat<npint>& arg) {return std::make_tuple(arg, arg);});
        
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif        
}

