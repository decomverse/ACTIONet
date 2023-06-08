#pragma once

#include <numpy/npy_common.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace std {

template <typename Num>
constexpr auto cbegin(const arma::Mat<Num>& matrix)
    -> decltype(std::begin(matrix)) {
  return std::begin(matrix);
}

template <typename Num>
constexpr auto cend(const arma::Mat<Num>& matrix)
    -> decltype(std::end(matrix)) {
  return std::end(matrix);
}

}  // namespace std

namespace aw {

using namespace std::literals::string_literals;
namespace py = pybind11;

using npulong = npy_ulong;
using arma::Col;
using arma::Cube;
using arma::Mat;
template <typename T>
using py_arr = py::array_t<T, py::array::f_style | py::array::forcecast>;
template <typename T>
using py_arr_c = py::array_t<T, py::array::c_style | py::array::forcecast>;
using npint = npy_int;
using npdouble = npy_double;
using intcube = arma::Cube<npint>;
using intmat = arma::Mat<npint>;
using dcube = arma::Cube<npdouble>;
using dmat = arma::Mat<npdouble>;
using dvec = arma::Col<npdouble>;
using iuvec = arma::Col<npulong>;
using std::begin;
using std::cbegin;
using std::cend;
using std::end;

/* Tests if U is a specialization of T */
template <template <typename...> typename T, typename U>
struct is_specialization_of : std::false_type {};

template <template <typename...> typename T, typename... Args>
struct is_specialization_of<T, T<Args...>> : std::true_type {};

template <typename T>
using is_SpMat = aw::is_specialization_of<arma::SpMat, T>;
template <typename T>
using is_Col = aw::is_specialization_of<arma::Col, T>;
template <typename T>
using is_Row = aw::is_specialization_of<arma::Row, T>;
template <typename T>
using is_Mat = aw::is_specialization_of<arma::Mat, T>;
template <typename T>
using is_Cube = aw::is_specialization_of<arma::Cube, T>;
template <typename T>
using is_Py_Arr = aw::is_specialization_of<py_arr, T>;

template <typename T1>
std::string stringify(const T1& arg1) {
  std::stringstream stream;
  stream << arg1;

  return stream.str();
}

template <typename T1, typename... Args>
std::string stringify(const T1& arg1, const Args&... args) {
  return stringify(arg1) + "\n"s + stringify(args...);
}

std::string stringify() { return ""s; }

template <typename Scalar, npulong dim, typename InfoType>
std::vector<npulong> extract_shape(const InfoType& info) {
  std::vector<npulong> shape(dim, 1);
  std::string format_err_msg =
      "The format descriptor strings are not the same. Are you using the right "
      " template specialization?";

  if (info.format != py::format_descriptor<Scalar>::format()) {
    // Temporary fix for sparse matrix
    if (info.format != "L" && py::format_descriptor<Scalar>::format() != "Q") {
      throw std::runtime_error(format_err_msg);
    }
  }

  if (sizeof(Scalar) != info.itemsize) {
    std::string size_err_msg =
        "The type you are storing the data in does not contain the same number "
        "of "
        "of bytes as the type you are storing the data in.";
    throw std::runtime_error(size_err_msg);
  }

  if (static_cast<npulong>(info.ndim) > dim) {
    throw std::runtime_error("Incompatible buffer dimensions");
  }
  std::copy(std::begin(info.shape), std::end(info.shape), std::begin(shape));

  return shape;
}

template <typename Scalar>
py::buffer_info def_buffer(arma::Row<Scalar>&& m) {
  py::buffer_info return_buf;

  if (m.size() >= 0) {
    return_buf = py::buffer_info(
        reinterpret_cast<void*>(m.memptr()), /* data */
        static_cast<ssize_t>(sizeof(Scalar)),
        py::format_descriptor<Scalar>::format(), static_cast<ssize_t>(1),
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_cols)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)) *
                             static_cast<ssize_t>(m.n_rows)});

  } else {
    return_buf = py::buffer_info(
        nullptr, /* data */
        static_cast<ssize_t>(sizeof(Scalar)),
        py::format_descriptor<Scalar>::format(), static_cast<ssize_t>(1),
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_cols)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)) *
                             static_cast<ssize_t>(m.n_rows)});
  }

  return return_buf;
}

template <typename Scalar>
py::buffer_info def_buffer(Col<Scalar>&& m) {
  if (m.size() > 0) {
    return py::buffer_info(
        reinterpret_cast<void*>(m.memptr()), /* data */
        static_cast<ssize_t>(sizeof(Scalar)),
        py::format_descriptor<Scalar>::format(), static_cast<ssize_t>(1),
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_rows)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar))});
  } else {
    return py::buffer_info(
        nullptr, static_cast<ssize_t>(sizeof(Scalar)),
        py::format_descriptor<Scalar>::format(), static_cast<ssize_t>(1),
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_rows)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar))});
  }
}

template <typename Scalar>
py::buffer_info def_buffer(Mat<Scalar>&& m) {
  py::buffer_info return_buf;

  if (m.size() > 0) {
    return_buf = py::buffer_info(
        reinterpret_cast<void*>(m.memptr()),     /* data */
        static_cast<ssize_t>(sizeof(Scalar)),    /* Size of one scalar. */
        py::format_descriptor<Scalar>::format(), /*Python struct-style format
                                                    descriptor */
        static_cast<ssize_t>(2),                 /*Number of dimensions */
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_rows),
                             static_cast<ssize_t>(m.n_cols)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)),
                             static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows)});

  } else {
    return_buf = py::buffer_info(
        nullptr,                                 /* Pointer to memory buffer. */
        static_cast<ssize_t>(sizeof(Scalar)),    /* Size of one scalar. */
        py::format_descriptor<Scalar>::format(), /*Python struct-style format
                                                    descriptor */
        static_cast<ssize_t>(2),                 /*Number of dimensions */
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_rows),
                             static_cast<ssize_t>(m.n_cols)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)),
                             static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows)});
  }

  return return_buf;
}

template <typename Scalar>
py::buffer_info def_buffer(Cube<Scalar>&& m) {
  py::buffer_info return_buf;

  if (m.size() >= 0) {
    return_buf = py::buffer_info(
        reinterpret_cast<void*>(m.memptr()),     /* data */
        static_cast<ssize_t>(sizeof(Scalar)),    /* Size of one scalar. */
        py::format_descriptor<Scalar>::format(), /*Python struct-style format
                                                    descriptor */
        static_cast<ssize_t>(3), {m.n_slices, m.n_rows, m.n_cols},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows) *
                                 static_cast<ssize_t>(m.n_cols),
                             static_cast<ssize_t>(sizeof(Scalar)),
                             static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows)});
  } else {
    return_buf = py::buffer_info(
        nullptr, static_cast<ssize_t>(sizeof(Scalar)), /* Size of one scalar. */
        py::format_descriptor<Scalar>::format(), /*Python struct-style format
                                                    descriptor */
        static_cast<ssize_t>(3),
        std::vector<ssize_t>{static_cast<ssize_t>(m.n_slices),
                             static_cast<ssize_t>(m.n_rows),
                             static_cast<ssize_t>(m.n_cols)},
        std::vector<ssize_t>{static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows) *
                                 static_cast<ssize_t>(m.n_cols),
                             static_cast<ssize_t>(sizeof(Scalar)),
                             static_cast<ssize_t>(sizeof(Scalar)) *
                                 static_cast<ssize_t>(m.n_rows)});
  }

  return return_buf;
}

/*
 * Converts something that may be an lvalue into an rvalue by copying if
 * necessary.
 */
template <typename T>
T clone(T val) {
  return val;
}

template <typename Scalar>
py::buffer_info def_buffer(const Mat<Scalar>& m) {
  return def_buffer<Scalar>(clone(m));
}

template <typename Scalar>
py::buffer_info def_buffer(const Cube<Scalar>& m) {
  return def_buffer<Scalar>(clone(m));
}

/*
 * These two functions call make_py_arr function which takes an Rvalue, and so
 * if they are given something else they copy it.
 */
template <template <typename> typename Arr, typename Scalar>
auto make_py_arr(Arr<Scalar>&& arr) {
  return py_arr<Scalar>(def_buffer(std::forward<Arr<Scalar>>(arr)));
}

template <template <typename> typename Arr, typename Scalar>
auto make_py_arr(const Arr<Scalar>& arr) {
  return py_arr<Scalar>(def_buffer(clone(arr)));
}

/*
 * This is a templated functor that has overloads that convert the various types
 * that I want to pass from Python to C++.
 */
template <typename ReturnType, typename SFINAE = std::true_type>
struct conv_to {
  static_assert(!SFINAE::value, "The general case is not defined.");
  template <typename InnerArgType>
  static ReturnType from(InnerArgType&&);
};

template <typename ReturnType>
struct conv_to<ReturnType, typename is_Mat<ReturnType>::type> {
  template <typename Scalar, int ForwardType>
  static ReturnType from(py::array_t<Scalar, ForwardType>& b) {
    if (ForwardType != py::array::f_style) {
      b = py::array_t<Scalar, py::array::f_style>(b);
    }

    py::buffer_info info;
    if (b.writeable()) {
      info = b.request(true);
    } else {
      info = b.request(false);
    }

    auto shape = extract_shape<Scalar, 2>(info);
    size_t length = shape[0] * shape[1];

    if (length > 0) {
      Scalar* memptr = static_cast<Scalar*>(info.ptr);

      if (b.writeable()) {
        return Mat<Scalar>(memptr, shape[0], shape[1]);

      } else {
        arma::Mat<Scalar> return_mat(shape[0], shape[1]);
        std::copy(memptr, memptr + length, std::begin(return_mat));

        return return_mat;
      }

    } else {
      return arma::Mat<Scalar>(shape[0], shape[1]);
    }
  }
};

/*
 * This function converts a 3 dimensional numpy array into an armadillo cube.
 * The tricky part is that the memory layout of the two arrays is different.
 * Consequently, you need to adjust the various paramters.
 */
template <typename ReturnType>
struct conv_to<ReturnType, typename is_Cube<ReturnType>::type> {
  template <typename Scalar, int ForwardType>
  static ReturnType from(py::array_t<Scalar, ForwardType>& b) {
    if (ForwardType != py::array::c_style) {
      b = py::array_t<Scalar, py::array::c_style>(b);
    }

    /* I copy the data and so I can get a non-writeable buffer */
    const auto info = b.request(false);
    const auto shape = extract_shape<Scalar, 3>(info);
    npulong length = std::accumulate(std::cbegin(shape), std::cend(shape), 1,
                                     std::multiplies<Scalar>());

    if (length > 0) {
      Cube<Scalar> tmp_mat(shape[2], shape[1], shape[0]);
      Cube<Scalar> return_mat(shape[1], shape[2], shape[0]);
      Scalar const* memptr = static_cast<Scalar*>(info.ptr);

      std::copy(memptr, memptr + length, std::begin(tmp_mat));

      for (npulong n = 0; n < return_mat.n_slices; ++n) {
        return_mat.slice(n) = tmp_mat.slice(n).t();
      }

      return return_mat;

    } else {
      return arma::Cube<Scalar>(shape[1], shape[2], shape[0]);
    }
  }
};

template <typename ReturnType>
struct conv_to<ReturnType, typename is_Row<ReturnType>::type> {
  template <typename Scalar, int ForwardType>
  static ReturnType from(py::array_t<Scalar, ForwardType>& b) {
    if (ForwardType != py::array::f_style &&
        ForwardType != py::array::c_style) {
      b = py::array_t<Scalar, py::array::c_style>(b);
    }

    if (b.writeable()) {
      auto info = b.request(true);
      auto shape = extract_shape<Scalar, 1>(info);

      if (shape[0] > 0) {
        return arma::Row<Scalar>(static_cast<Scalar*>(info.ptr), shape[0]);

      } else {
        return arma::Row<Scalar>(shape[0]);
      }
    } else {
      auto info = b.request(false);
      auto shape = extract_shape<Scalar, 1>(info);
      auto ref_ptr = static_cast<Scalar*>(info.ptr);
      arma::Row<Scalar> return_row(shape[0]);
      std::copy(ref_ptr, ref_ptr + shape[0], std::begin(return_row));
      return return_row;
    }
  }
};

template <typename ReturnType>
struct conv_to<ReturnType, typename is_Col<ReturnType>::type> {
  template <typename Scalar, int ForwardType>
  static ReturnType from(py::array_t<Scalar, ForwardType>& b) {
    if (ForwardType != py::array::f_style &&
        ForwardType != py::array::c_style) {
      b = py::array_t<Scalar, py::array::c_style>(b);
    }

    if (b.writeable()) {
      auto info = b.request(true);
      auto shape = extract_shape<Scalar, 1>(info);

      if (shape[0] > 0) {
        return arma::Col<Scalar>(static_cast<Scalar*>(info.ptr), shape[0]);

      } else {
        return arma::Col<Scalar>(shape[0]);
      }

    } else {
      auto info = b.request(false);
      auto shape = extract_shape<Scalar, 1>(info);

      arma::Col<Scalar> return_mat(shape[0]);
      auto ref_ptr = static_cast<Scalar*>(info.ptr);
      std::copy(ref_ptr, ref_ptr + shape[0], std::begin(return_mat));
      return return_mat;
    }
  }
};

template <typename T1, typename T2>
bool check_same_n_rows(const T1& first, const T2& second) {
  return (first.n_rows == second.n_rows);
}

template <typename T1, typename T2, typename... Args>
bool check_same_n_rows(const T1& first, const T2& second, const Args&... args) {
  return check_same_n_rows(first, second) && check_same_n_rows(second, args...);
}

template <typename... Args>
void assert_same_n_rows(const Args&... args) {
  assert(check_same_n_rows(args...));
}

template <typename... Args>
void claim_same_n_rows(const Args&... args) {
  if (!check_same_n_rows(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(
        "The number of rows in each of the "s + std::to_string(num_args) +
        " arrays is not the same."s + " The args are"s + stringify(args...));
  }
}

template <typename T1, typename T2>
bool check_same_n_cols(const T1& first, const T2& second) {
  return (first.n_cols == second.n_cols);
}

template <typename T1, typename T2, typename... Args>
bool check_same_n_cols(const T1& first, const T2& second, const Args&... args) {
  return check_same_n_cols(first, second) && check_same_n_cols(second, args...);
}

template <typename... Args>
void assert_same_n_cols(const Args&... args) {
  assert(check_same_n_cols(args...));
}

template <typename... Args>
void claim_same_n_cols(const Args&... args) {
  if (!check_same_n_cols(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(
        "The number of cols in each of the "s + std::to_string(num_args) +
        " arrays is not the same."s + " The args are"s + stringify(args...));
  }
}

template <typename T1, typename T2>
bool check_same_n_slices(const T1& first, const T2& second) {
  return (first.n_slices == second.n_slices);
}

template <typename T1, typename T2, typename... Args>
bool check_same_n_slices(const T1& first, const T2& second,
                         const Args&... args) {
  return check_same_n_slices(first, second) &&
         check_same_n_slices(second, args...);
}

template <typename... Args>
void assert_same_n_slices(const Args&... args) {
  assert(check_same_n_slices(args...));
}

template <typename... Args>
void claim_same_n_slices(const Args&... args) {
  if (!check_same_n_slices(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(
        "The number of slices in each of the "s + std::to_string(num_args) +
        " arrays is not the same."s + " The args are"s + stringify(args...));
  }
}

template <typename T1, typename T2>
bool check_same_size(const T1& first, const T2& second) {
  return (arma::size(first) == arma::size(second));
}

template <typename T1, typename T2, typename... Args>
bool check_same_size(const T1& first, const T2& second, const Args&... args) {
  return check_same_size(first, second) && check_same_size(second, args...);
}

template <typename... Args>
void assert_same_size(const Args&... args) {
  assert(check_same_size(args...));
}

template <typename... Args>
void claim_same_size(const Args&... args) {
  if (!check_same_size(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(
        std::string("All of the ") + std::to_string(num_args) +
        std::string(" arrays do not have the same shape."));
  }
}

template <typename T1, typename T2>
bool check_same(const T1& first, const T2& second) {
  return (first == second);
}

template <typename T1, typename T2, typename... Args>
bool check_same(const T1& first, const T2& second, const Args&... args) {
  return check_same(first, second) && check_same(second, args...);
}

template <typename... Args>
void assert_same(const Args&... args) {
  assert(check_same(args...));
}

template <typename... Args>
void claim_same(const Args&... args) {
  if (!check_same(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(
        std::string("All of the ") + std::to_string(num_args) +
        std::string(" arguments do not equal each other."));
  }
}

template <typename T1>
bool check_is_finite(const T1& first) {
  return first.is_finite();
}

template <typename T1, typename... Args>
bool check_is_finite(const T1& first, const Args&... args) {
  return check_is_finite(first) && check_is_finite(args...);
}

template <typename... Args>
void assert_is_finite(const Args&... args) {
  assert(check_is_finite(args...));
}

template <typename... Args>
void claim_is_finite(const Args&... args) {
  if (!check_is_finite(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(std::string("Not all of the ") +
                             std::to_string(num_args) +
                             std::string(" arrays are entirely finite."));
  }
}

template <typename T1>
bool check_square(const T1& first) {
  return (first.n_cols == first.n_rows);
}

template <typename T1, typename... Args>
bool check_square(const T1& first, const Args&... args) {
  return check_square(first) && check_square(args...);
}

template <typename... Args>
void assert_square(const Args&... args) {
  assert(check_square(args...));
}

template <typename... Args>
void claim_square(const Args&... args) {
  if (!check_squre(args...)) {
    std::size_t num_args = sizeof...(args);
    throw std::runtime_error(std::string("Not all of the ") +
                             std::to_string(num_args) +
                             std::string(" arrays are square."));
  }
}

/*
 * Takes a vector and appends an element to it.
 */
template <typename NumericType>
Col<NumericType> append(Col<NumericType> arr, const NumericType val) {
  arr.resize(arr.n_elem + 1);
  arr.tail(1) = val;

  return arr;
}

template <typename NumericType>
arma::Row<NumericType> append(arma::Row<NumericType> arr,
                              const NumericType val) {
  arr.resize(arr.n_elem + 1);
  arr.tail(1) = val;

  return arr;
}

/* Takes a matrix and appends a row to it. */
template <typename NumericType>
Mat<NumericType> append(Mat<NumericType> arr,
                        const arma::Row<NumericType> row) {
  assert_same_n_cols(arr, row);
  arr.resize(arr.n_rows + 1, arr.n_cols);
  arr.tail_rows(1) = row;

  return arr;
}

/*Takes a matrix and appends a column to it. */
template <typename NumericType>
arma::Mat<NumericType> append(Mat<NumericType> arr,
                              const Col<NumericType> vec) {
  assert_same_n_rows(arr, vec);
  arr.resize(arr.n_rows, arr.n_cols + 1);
  arr.tail_cols(1) = vec;

  return arr;
}

/* Appends a matrix to the end of a cube. */
template <typename NumericType>
arma::Cube<NumericType> append(arma::Cube<NumericType> arr,
                               const Mat<NumericType> mat) {
  assert_same_n_cols(arr, mat);
  assert_same_n_rows(arr, mat);

  arr.resize(arr.n_rows, arr.n_cols, arr.n_slices + 1);
  arr.tail_slice(1) = mat;
}

template <typename NumericType>
void append_inplace(Col<NumericType>& arr, const NumericType val) {
  arr.resize(arr.n_elem + 1);
  arr.tail(1) = val;
}

template <typename NumericType>
void append_inplace(arma::Row<NumericType>& arr, const NumericType val) {
  arr.resize(arr.n_elem + 1);
  arr.tail(1) = val;
}

/* Takes a matrix and appends a row to it. */
template <typename NumericType>
void append_inplace(Mat<NumericType>& arr, const arma::Row<NumericType> row) {
  assert_same_n_cols(arr, row);
  arr.resize(arr.n_rows + 1, arr.n_cols);
  arr.tail_rows(1) = row;
}

/*Takes a matrix and appends a column to it. */
template <typename NumericType>
void append_inplace(Mat<NumericType>& arr, const Col<NumericType> vec) {
  assert_same_n_rows(arr, vec);
  arr.resize(arr.n_rows, arr.n_cols + 1);
  arr.tail_cols(1) = vec;
}

/* Appends a matrix to the end of a cube. */
template <typename NumericType>
void append_inplace(arma::Cube<NumericType>& arr, const Mat<NumericType> mat) {
  assert_same_n_cols(arr, mat);
  assert_same_n_rows(arr, mat);

  arr.resize(arr.n_rows, arr.n_cols, arr.n_slices + 1);
  arr.tail_slices(1) = mat;
}

/*
 * Converts a container to a scalar. It calls the as_scalar funciton in
 * Armadillo if it exists.
 */
template <typename ContainerType>
auto as_scalar(ContainerType&& container)
    -> decltype(arma::as_scalar(std::forward<ContainerType>(container))) {
  return arma::as_scalar(std::forward<ContainerType>(container));
}

/*
 * I provide a template specialization for std::vector. I convert to an
 * Armadillo vector and then call the std::vector there. By doing this I can
 * recover the same exceptions as in the case above.
 */
template <typename ValueType>
ValueType as_scalar(const std::vector<ValueType>& container) {
  return arma::as_scalar(Col<ValueType>(container));
}

/*
 * Converts a vector to its square matrix as long as the vector has a Scalarber
 * of elmeents that is a perfect square.
 */
template <typename ValueType>
Mat<ValueType> reshape_square(Col<ValueType>&& vec) {
  npulong side_length = (npulong)std::sqrt(vec.n_elem);
  assert(side_length * side_length == vec.n_elem &&
         "The vector does not have a perfect square Scalarber of "
         "elements");
  return arma::reshape(std::forward<Col<ValueType>>(vec), side_length,
                       side_length);
}

template <typename ValueType>
Mat<ValueType> reshape_square(const std::vector<ValueType>& vec) {
  return reshape_square(Col<ValueType>(vec));
}

}  // namespace aw

namespace pybind11 {
namespace detail {

using std::enable_if_t;

/*
 * Handles conversions for dense Armadillo arrays.
 */
template <typename Type>
struct type_caster<
    Type, enable_if_t<(aw::is_Col<Type>::value || aw::is_Row<Type>::value) ||
                      (aw::is_Mat<Type>::value || aw::is_Cube<Type>::value)>> {
  using Scalar = typename std::decay_t<typename Type::elem_type>;

  /*
   * This function loads python arrays and coverts them to Armadillo columns.
   */

  bool load(handle src, bool convert) {
    if (!convert && !isinstance<array_t<Scalar>>(src)) {
      return false;
    }

    auto buf = array_t<Scalar>::ensure(src);

    if (!buf) {
      return false;
    }

    value = aw::conv_to<Type>::from(buf);

    return true;
  }

  static handle cast(const Type& src, return_value_policy /* policy */,
                     handle /* parent */) {
    auto returnval = aw::make_py_arr(src).release();

    return returnval;
  }

  PYBIND11_TYPE_CASTER(Type, _("Numpy.ndarray[") +
                                 npy_format_descriptor<Scalar>::name + _("]"));
};

template <typename Type>
struct type_caster<Type, enable_if_t<aw::is_SpMat<Type>::value>> {
  using Scalar = typename std::decay_t<typename Type::elem_type>;
  using StorageIndexType =
      typename std::decay_t<decltype(*std::declval<Type>().row_indices)>;
  using IndexType =
      typename std::decay_t<decltype(std::declval<Type>().n_nonzero)>;

  using StorageIndexVector = arma::Col<StorageIndexType>;
  using IndvecVector = arma::Col<IndexType>;
  using DataMat = arma::Mat<Scalar>;

  /*
   * This function transfers python sparse arrays to armadillo. Armadillo only
   * accepts csc format matrices, and do it converts them if necessary.
   */
  bool load(handle src, bool) {
    if (!src) {
      return false;
    }

    handle obj = reinterpret_borrow<object>(src);
    object sparse_module = module::import("scipy.sparse");
    object matrix_type = sparse_module.attr("csc_matrix");

    if (!obj.get_type().is(matrix_type)) {
      try {
        obj = matrix_type(obj);
      } catch (const error_already_set&) {
        return false;
      }
    }

    auto values = array_t<Scalar>((object)obj.attr("data"));
    auto row_indices = array_t<StorageIndexType>((object)obj.attr("indices"));
    auto col_ptrs = array_t<StorageIndexType>((object)obj.attr("indptr"));
    auto shape = pybind11::tuple((pybind11::object)obj.attr("shape"));

    if (!values || !row_indices || !col_ptrs) {
      return false;
    }

    arma::uvec row_vector = arma::conv_to<arma::uvec>::from(
        aw::conv_to<StorageIndexVector>::from(row_indices));
    arma::uvec col_ptr_vector = arma::conv_to<arma::uvec>::from(
        aw::conv_to<StorageIndexVector>::from(col_ptrs));

    arma::Mat<Scalar> dmat_values = aw::conv_to<DataMat>::from(values);
    arma::Col<Scalar> values_vector = arma::vectorise(dmat_values);
    value = arma::SpMat<Scalar>(row_vector, col_ptr_vector, values_vector,
                                shape[0].cast<IndexType>(),
                                shape[1].cast<IndexType>());

    return true;
  }

  /*
   * This function transfers the data from the csc format in armadillo to the
   * equivalenet format in python.
   */
  static handle cast(const Type& src, return_value_policy /* policy */,
                     handle /* parent */) {
    object matrix_type = module::import("scipy.sparse").attr("csc_matrix");

    // Armadillo pads the arrays with an extra elements, and so we only copy the
    // elements that are necessary in Python.
    array data(static_cast<size_t>(src.n_nonzero), src.values);
    array col_ptrs(static_cast<size_t>(src.n_cols + 1), src.col_ptrs);
    array row_indices(static_cast<size_t>(src.n_nonzero), src.row_indices);

    return matrix_type(std::make_tuple(data, row_indices, col_ptrs),
                       std::make_pair(static_cast<aw::npulong>(src.n_rows),
                                      static_cast<aw::npulong>(src.n_cols)))
        .release();
  }

  PYBIND11_TYPE_CASTER(Type,
                       _("scipy.sparse.csc_matrix[") +
                           npy_format_descriptor<std::decay_t<Scalar>>::name +
                           _("]"));
 };
}
}
