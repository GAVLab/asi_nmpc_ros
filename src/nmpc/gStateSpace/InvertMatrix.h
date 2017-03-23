
/**************************************************************************************
 *
 *  InvertMatrix.h
 *
 *  Mike Newey
 *  May 2006
 *
 *  This code is from the ublas website (freely available) and will invert a ublas matrix.
 *				 
 *
 *********************************************************************************************/



 #ifndef _INVERT_MATRIX_H_
 #define _INVERT_MATRIX_H_

 #undef BOOST_UBLAS_TYPE_CHECK
 #define BOOST_UBLAS_TYPE_CHECK 0


 // REMEMBER to update "lu.hpp" header includes from boost-CVS
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/vector_proxy.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/triangular.hpp>
 #include <boost/numeric/ublas/lu.hpp>
 #include <boost/numeric/ublas/io.hpp>
 #include <stdio.h>

 

 namespace ublas = boost::numeric::ublas;

 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 template<class T>
 int InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
        //check to make sure the matrix is square and invertible
        if (input.size1() != input.size2())
            return -21;
               
 	// create a working copy of the input
 	matrix<T> A(input);
        // matrix<T> A;
        // A.assign(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());

        // std::cout << "invertmatrix" << std::endl;
        //std::cout << A << std::endl;
 	// perform LU-factorization
	int res = lu_factorize(A, pm);	
 	if( res != 0 ) return -22;

 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);

        return 0;
 }

 #undef BOOST_UBLAS_TYPE_CHECK
 #define BOOST_UBLAS_TYPE_CHECK
 #endif //_INVERT_MATRIX_H_


//Error Codes:
// -1:  Not a square matrix
// -2:  Singular Matrix


