#ifndef GSSCOMMON_H
#define GSSCOMMON_H

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<double> UDvector;
typedef boost::numeric::ublas::matrix<double> UDmatrix;
typedef boost::numeric::ublas::vector<int> UIvector;
typedef boost::numeric::ublas::matrix<int> UImatrix;
typedef boost::numeric::ublas::identity_matrix<double> IDmatrix;
typedef boost::numeric::ublas::triangular_matrix<double, boost::numeric::ublas::upper, boost::numeric::ublas::row_major> UTriMatrix;

#include "InvertMatrix.h"

#endif