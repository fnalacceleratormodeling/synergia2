var data = {lines:[
{"lineNum":"    1","line":"// This file is part of Eigen, a lightweight C++ template library"},
{"lineNum":"    2","line":"// for linear algebra."},
{"lineNum":"    3","line":"//"},
{"lineNum":"    4","line":"// Copyright (C) 2008-2010 Benoit Jacob <jacob.benoit.1@gmail.com>"},
{"lineNum":"    5","line":"// Copyright (C) 2014 Gael Guennebaud <gael.guennebaud@inria.fr>"},
{"lineNum":"    6","line":"//"},
{"lineNum":"    7","line":"// This Source Code Form is subject to the terms of the Mozilla"},
{"lineNum":"    8","line":"// Public License v. 2.0. If a copy of the MPL was not distributed"},
{"lineNum":"    9","line":"// with this file, You can obtain one at http://mozilla.org/MPL/2.0/."},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#ifndef EIGEN_INVERSE_IMPL_H"},
{"lineNum":"   12","line":"#define EIGEN_INVERSE_IMPL_H"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#include \"./InternalHeaderCheck.h\""},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"namespace Eigen {"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"namespace internal {"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"/**********************************"},
{"lineNum":"   21","line":"*** General case implementation ***"},
{"lineNum":"   22","line":"**********************************/"},
{"lineNum":"   23","line":""},
{"lineNum":"   24","line":"template<typename MatrixType, typename ResultType, int Size = MatrixType::RowsAtCompileTime>"},
{"lineNum":"   25","line":"struct compute_inverse"},
{"lineNum":"   26","line":"{"},
{"lineNum":"   27","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"   28","line":"  static inline void run(const MatrixType& matrix, ResultType& result)"},
{"lineNum":"   29","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   30","line":"    result = matrix.partialPivLu().inverse();","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   31","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"   32","line":"};"},
{"lineNum":"   33","line":""},
{"lineNum":"   34","line":"template<typename MatrixType, typename ResultType, int Size = MatrixType::RowsAtCompileTime>"},
{"lineNum":"   35","line":"struct compute_inverse_and_det_with_check { /* nothing! general case not supported. */ };"},
{"lineNum":"   36","line":""},
{"lineNum":"   37","line":"/****************************"},
{"lineNum":"   38","line":"*** Size 1 implementation ***"},
{"lineNum":"   39","line":"****************************/"},
{"lineNum":"   40","line":""},
{"lineNum":"   41","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"   42","line":"struct compute_inverse<MatrixType, ResultType, 1>"},
{"lineNum":"   43","line":"{"},
{"lineNum":"   44","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"   45","line":"  static inline void run(const MatrixType& matrix, ResultType& result)"},
{"lineNum":"   46","line":"  {"},
{"lineNum":"   47","line":"    typedef typename MatrixType::Scalar Scalar;"},
{"lineNum":"   48","line":"    internal::evaluator<MatrixType> matrixEval(matrix);"},
{"lineNum":"   49","line":"    result.coeffRef(0,0) = Scalar(1) / matrixEval.coeff(0,0);"},
{"lineNum":"   50","line":"  }"},
{"lineNum":"   51","line":"};"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"   54","line":"struct compute_inverse_and_det_with_check<MatrixType, ResultType, 1>"},
{"lineNum":"   55","line":"{"},
{"lineNum":"   56","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"   57","line":"  static inline void run("},
{"lineNum":"   58","line":"    const MatrixType& matrix,"},
{"lineNum":"   59","line":"    const typename MatrixType::RealScalar& absDeterminantThreshold,"},
{"lineNum":"   60","line":"    ResultType& result,"},
{"lineNum":"   61","line":"    typename ResultType::Scalar& determinant,"},
{"lineNum":"   62","line":"    bool& invertible"},
{"lineNum":"   63","line":"  )"},
{"lineNum":"   64","line":"  {"},
{"lineNum":"   65","line":"    using std::abs;"},
{"lineNum":"   66","line":"    determinant = matrix.coeff(0,0);"},
{"lineNum":"   67","line":"    invertible = abs(determinant) > absDeterminantThreshold;"},
{"lineNum":"   68","line":"    if(invertible) result.coeffRef(0,0) = typename ResultType::Scalar(1) / determinant;"},
{"lineNum":"   69","line":"  }"},
{"lineNum":"   70","line":"};"},
{"lineNum":"   71","line":""},
{"lineNum":"   72","line":"/****************************"},
{"lineNum":"   73","line":"*** Size 2 implementation ***"},
{"lineNum":"   74","line":"****************************/"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"   77","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"   78","line":"inline void compute_inverse_size2_helper("},
{"lineNum":"   79","line":"    const MatrixType& matrix, const typename ResultType::Scalar& invdet,"},
{"lineNum":"   80","line":"    ResultType& result)"},
{"lineNum":"   81","line":"{"},
{"lineNum":"   82","line":"  typename ResultType::Scalar temp = matrix.coeff(0,0);"},
{"lineNum":"   83","line":"  result.coeffRef(0,0) =  matrix.coeff(1,1) * invdet;"},
{"lineNum":"   84","line":"  result.coeffRef(1,0) = -matrix.coeff(1,0) * invdet;"},
{"lineNum":"   85","line":"  result.coeffRef(0,1) = -matrix.coeff(0,1) * invdet;"},
{"lineNum":"   86","line":"  result.coeffRef(1,1) =  temp * invdet;"},
{"lineNum":"   87","line":"}"},
{"lineNum":"   88","line":""},
{"lineNum":"   89","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"   90","line":"struct compute_inverse<MatrixType, ResultType, 2>"},
{"lineNum":"   91","line":"{"},
{"lineNum":"   92","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"   93","line":"  static inline void run(const MatrixType& matrix, ResultType& result)"},
{"lineNum":"   94","line":"  {"},
{"lineNum":"   95","line":"    typedef typename ResultType::Scalar Scalar;"},
{"lineNum":"   96","line":"    const Scalar invdet = typename MatrixType::Scalar(1) / matrix.determinant();"},
{"lineNum":"   97","line":"    compute_inverse_size2_helper(matrix, invdet, result);"},
{"lineNum":"   98","line":"  }"},
{"lineNum":"   99","line":"};"},
{"lineNum":"  100","line":""},
{"lineNum":"  101","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  102","line":"struct compute_inverse_and_det_with_check<MatrixType, ResultType, 2>"},
{"lineNum":"  103","line":"{"},
{"lineNum":"  104","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  105","line":"  static inline void run("},
{"lineNum":"  106","line":"    const MatrixType& matrix,"},
{"lineNum":"  107","line":"    const typename MatrixType::RealScalar& absDeterminantThreshold,"},
{"lineNum":"  108","line":"    ResultType& inverse,"},
{"lineNum":"  109","line":"    typename ResultType::Scalar& determinant,"},
{"lineNum":"  110","line":"    bool& invertible"},
{"lineNum":"  111","line":"  )"},
{"lineNum":"  112","line":"  {"},
{"lineNum":"  113","line":"    using std::abs;"},
{"lineNum":"  114","line":"    typedef typename ResultType::Scalar Scalar;"},
{"lineNum":"  115","line":"    determinant = matrix.determinant();"},
{"lineNum":"  116","line":"    invertible = abs(determinant) > absDeterminantThreshold;"},
{"lineNum":"  117","line":"    if(!invertible) return;"},
{"lineNum":"  118","line":"    const Scalar invdet = Scalar(1) / determinant;"},
{"lineNum":"  119","line":"    compute_inverse_size2_helper(matrix, invdet, inverse);"},
{"lineNum":"  120","line":"  }"},
{"lineNum":"  121","line":"};"},
{"lineNum":"  122","line":""},
{"lineNum":"  123","line":"/****************************"},
{"lineNum":"  124","line":"*** Size 3 implementation ***"},
{"lineNum":"  125","line":"****************************/"},
{"lineNum":"  126","line":""},
{"lineNum":"  127","line":"template<typename MatrixType, int i, int j>"},
{"lineNum":"  128","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"  129","line":"inline typename MatrixType::Scalar cofactor_3x3(const MatrixType& m)"},
{"lineNum":"  130","line":"{"},
{"lineNum":"  131","line":"  enum {"},
{"lineNum":"  132","line":"    i1 = (i+1) % 3,"},
{"lineNum":"  133","line":"    i2 = (i+2) % 3,"},
{"lineNum":"  134","line":"    j1 = (j+1) % 3,"},
{"lineNum":"  135","line":"    j2 = (j+2) % 3"},
{"lineNum":"  136","line":"  };"},
{"lineNum":"  137","line":"  return m.coeff(i1, j1) * m.coeff(i2, j2)"},
{"lineNum":"  138","line":"       - m.coeff(i1, j2) * m.coeff(i2, j1);"},
{"lineNum":"  139","line":"}"},
{"lineNum":"  140","line":""},
{"lineNum":"  141","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  142","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"  143","line":"inline void compute_inverse_size3_helper("},
{"lineNum":"  144","line":"    const MatrixType& matrix,"},
{"lineNum":"  145","line":"    const typename ResultType::Scalar& invdet,"},
{"lineNum":"  146","line":"    const Matrix<typename ResultType::Scalar,3,1>& cofactors_col0,"},
{"lineNum":"  147","line":"    ResultType& result)"},
{"lineNum":"  148","line":"{"},
{"lineNum":"  149","line":"  // Compute cofactors in a way that avoids aliasing issues."},
{"lineNum":"  150","line":"  typedef typename ResultType::Scalar Scalar;"},
{"lineNum":"  151","line":"  const Scalar c01 = cofactor_3x3<MatrixType,0,1>(matrix) * invdet;"},
{"lineNum":"  152","line":"  const Scalar c11 = cofactor_3x3<MatrixType,1,1>(matrix) * invdet;"},
{"lineNum":"  153","line":"  const Scalar c02 = cofactor_3x3<MatrixType,0,2>(matrix) * invdet;"},
{"lineNum":"  154","line":"  result.coeffRef(1,2) =  cofactor_3x3<MatrixType,2,1>(matrix) * invdet;"},
{"lineNum":"  155","line":"  result.coeffRef(2,1) =  cofactor_3x3<MatrixType,1,2>(matrix) * invdet;"},
{"lineNum":"  156","line":"  result.coeffRef(2,2) =  cofactor_3x3<MatrixType,2,2>(matrix) * invdet;"},
{"lineNum":"  157","line":"  result.coeffRef(1,0) =  c01;"},
{"lineNum":"  158","line":"  result.coeffRef(1,1) =  c11;"},
{"lineNum":"  159","line":"  result.coeffRef(2,0) =  c02;"},
{"lineNum":"  160","line":"  result.row(0) = cofactors_col0 * invdet;"},
{"lineNum":"  161","line":"}"},
{"lineNum":"  162","line":""},
{"lineNum":"  163","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  164","line":"struct compute_inverse<MatrixType, ResultType, 3>"},
{"lineNum":"  165","line":"{"},
{"lineNum":"  166","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  167","line":"  static inline void run(const MatrixType& matrix, ResultType& result)"},
{"lineNum":"  168","line":"  {"},
{"lineNum":"  169","line":"    typedef typename ResultType::Scalar Scalar;"},
{"lineNum":"  170","line":"    Matrix<typename MatrixType::Scalar,3,1> cofactors_col0;"},
{"lineNum":"  171","line":"    cofactors_col0.coeffRef(0) =  cofactor_3x3<MatrixType,0,0>(matrix);"},
{"lineNum":"  172","line":"    cofactors_col0.coeffRef(1) =  cofactor_3x3<MatrixType,1,0>(matrix);"},
{"lineNum":"  173","line":"    cofactors_col0.coeffRef(2) =  cofactor_3x3<MatrixType,2,0>(matrix);"},
{"lineNum":"  174","line":"    const Scalar det = (cofactors_col0.cwiseProduct(matrix.col(0))).sum();"},
{"lineNum":"  175","line":"    const Scalar invdet = Scalar(1) / det;"},
{"lineNum":"  176","line":"    compute_inverse_size3_helper(matrix, invdet, cofactors_col0, result);"},
{"lineNum":"  177","line":"  }"},
{"lineNum":"  178","line":"};"},
{"lineNum":"  179","line":""},
{"lineNum":"  180","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  181","line":"struct compute_inverse_and_det_with_check<MatrixType, ResultType, 3>"},
{"lineNum":"  182","line":"{"},
{"lineNum":"  183","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  184","line":"  static inline void run("},
{"lineNum":"  185","line":"    const MatrixType& matrix,"},
{"lineNum":"  186","line":"    const typename MatrixType::RealScalar& absDeterminantThreshold,"},
{"lineNum":"  187","line":"    ResultType& inverse,"},
{"lineNum":"  188","line":"    typename ResultType::Scalar& determinant,"},
{"lineNum":"  189","line":"    bool& invertible"},
{"lineNum":"  190","line":"  )"},
{"lineNum":"  191","line":"  {"},
{"lineNum":"  192","line":"    typedef typename ResultType::Scalar Scalar;"},
{"lineNum":"  193","line":"    Matrix<Scalar,3,1> cofactors_col0;"},
{"lineNum":"  194","line":"    cofactors_col0.coeffRef(0) =  cofactor_3x3<MatrixType,0,0>(matrix);"},
{"lineNum":"  195","line":"    cofactors_col0.coeffRef(1) =  cofactor_3x3<MatrixType,1,0>(matrix);"},
{"lineNum":"  196","line":"    cofactors_col0.coeffRef(2) =  cofactor_3x3<MatrixType,2,0>(matrix);"},
{"lineNum":"  197","line":"    determinant = (cofactors_col0.cwiseProduct(matrix.col(0))).sum();"},
{"lineNum":"  198","line":"    invertible = Eigen::numext::abs(determinant) > absDeterminantThreshold;"},
{"lineNum":"  199","line":"    if(!invertible) return;"},
{"lineNum":"  200","line":"    const Scalar invdet = Scalar(1) / determinant;"},
{"lineNum":"  201","line":"    compute_inverse_size3_helper(matrix, invdet, cofactors_col0, inverse);"},
{"lineNum":"  202","line":"  }"},
{"lineNum":"  203","line":"};"},
{"lineNum":"  204","line":""},
{"lineNum":"  205","line":"/****************************"},
{"lineNum":"  206","line":"*** Size 4 implementation ***"},
{"lineNum":"  207","line":"****************************/"},
{"lineNum":"  208","line":""},
{"lineNum":"  209","line":"template<typename Derived>"},
{"lineNum":"  210","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"  211","line":"inline const typename Derived::Scalar general_det3_helper"},
{"lineNum":"  212","line":"(const MatrixBase<Derived>& matrix, int i1, int i2, int i3, int j1, int j2, int j3)"},
{"lineNum":"  213","line":"{"},
{"lineNum":"  214","line":"  return matrix.coeff(i1,j1)"},
{"lineNum":"  215","line":"         * (matrix.coeff(i2,j2) * matrix.coeff(i3,j3) - matrix.coeff(i2,j3) * matrix.coeff(i3,j2));"},
{"lineNum":"  216","line":"}"},
{"lineNum":"  217","line":""},
{"lineNum":"  218","line":"template<typename MatrixType, int i, int j>"},
{"lineNum":"  219","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"  220","line":"inline typename MatrixType::Scalar cofactor_4x4(const MatrixType& matrix)"},
{"lineNum":"  221","line":"{"},
{"lineNum":"  222","line":"  enum {"},
{"lineNum":"  223","line":"    i1 = (i+1) % 4,"},
{"lineNum":"  224","line":"    i2 = (i+2) % 4,"},
{"lineNum":"  225","line":"    i3 = (i+3) % 4,"},
{"lineNum":"  226","line":"    j1 = (j+1) % 4,"},
{"lineNum":"  227","line":"    j2 = (j+2) % 4,"},
{"lineNum":"  228","line":"    j3 = (j+3) % 4"},
{"lineNum":"  229","line":"  };"},
{"lineNum":"  230","line":"  return general_det3_helper(matrix, i1, i2, i3, j1, j2, j3)"},
{"lineNum":"  231","line":"       + general_det3_helper(matrix, i2, i3, i1, j1, j2, j3)"},
{"lineNum":"  232","line":"       + general_det3_helper(matrix, i3, i1, i2, j1, j2, j3);"},
{"lineNum":"  233","line":"}"},
{"lineNum":"  234","line":""},
{"lineNum":"  235","line":"template<int Arch, typename Scalar, typename MatrixType, typename ResultType>"},
{"lineNum":"  236","line":"struct compute_inverse_size4"},
{"lineNum":"  237","line":"{"},
{"lineNum":"  238","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  239","line":"  static void run(const MatrixType& matrix, ResultType& result)"},
{"lineNum":"  240","line":"  {"},
{"lineNum":"  241","line":"    result.coeffRef(0,0) =  cofactor_4x4<MatrixType,0,0>(matrix);"},
{"lineNum":"  242","line":"    result.coeffRef(1,0) = -cofactor_4x4<MatrixType,0,1>(matrix);"},
{"lineNum":"  243","line":"    result.coeffRef(2,0) =  cofactor_4x4<MatrixType,0,2>(matrix);"},
{"lineNum":"  244","line":"    result.coeffRef(3,0) = -cofactor_4x4<MatrixType,0,3>(matrix);"},
{"lineNum":"  245","line":"    result.coeffRef(0,2) =  cofactor_4x4<MatrixType,2,0>(matrix);"},
{"lineNum":"  246","line":"    result.coeffRef(1,2) = -cofactor_4x4<MatrixType,2,1>(matrix);"},
{"lineNum":"  247","line":"    result.coeffRef(2,2) =  cofactor_4x4<MatrixType,2,2>(matrix);"},
{"lineNum":"  248","line":"    result.coeffRef(3,2) = -cofactor_4x4<MatrixType,2,3>(matrix);"},
{"lineNum":"  249","line":"    result.coeffRef(0,1) = -cofactor_4x4<MatrixType,1,0>(matrix);"},
{"lineNum":"  250","line":"    result.coeffRef(1,1) =  cofactor_4x4<MatrixType,1,1>(matrix);"},
{"lineNum":"  251","line":"    result.coeffRef(2,1) = -cofactor_4x4<MatrixType,1,2>(matrix);"},
{"lineNum":"  252","line":"    result.coeffRef(3,1) =  cofactor_4x4<MatrixType,1,3>(matrix);"},
{"lineNum":"  253","line":"    result.coeffRef(0,3) = -cofactor_4x4<MatrixType,3,0>(matrix);"},
{"lineNum":"  254","line":"    result.coeffRef(1,3) =  cofactor_4x4<MatrixType,3,1>(matrix);"},
{"lineNum":"  255","line":"    result.coeffRef(2,3) = -cofactor_4x4<MatrixType,3,2>(matrix);"},
{"lineNum":"  256","line":"    result.coeffRef(3,3) =  cofactor_4x4<MatrixType,3,3>(matrix);"},
{"lineNum":"  257","line":"    result /= (matrix.col(0).cwiseProduct(result.row(0).transpose())).sum();"},
{"lineNum":"  258","line":"  }"},
{"lineNum":"  259","line":"};"},
{"lineNum":"  260","line":""},
{"lineNum":"  261","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  262","line":"struct compute_inverse<MatrixType, ResultType, 4>"},
{"lineNum":"  263","line":" : compute_inverse_size4<Architecture::Target, typename MatrixType::Scalar,"},
{"lineNum":"  264","line":"                            MatrixType, ResultType>"},
{"lineNum":"  265","line":"{"},
{"lineNum":"  266","line":"};"},
{"lineNum":"  267","line":""},
{"lineNum":"  268","line":"template<typename MatrixType, typename ResultType>"},
{"lineNum":"  269","line":"struct compute_inverse_and_det_with_check<MatrixType, ResultType, 4>"},
{"lineNum":"  270","line":"{"},
{"lineNum":"  271","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  272","line":"  static inline void run("},
{"lineNum":"  273","line":"    const MatrixType& matrix,"},
{"lineNum":"  274","line":"    const typename MatrixType::RealScalar& absDeterminantThreshold,"},
{"lineNum":"  275","line":"    ResultType& inverse,"},
{"lineNum":"  276","line":"    typename ResultType::Scalar& determinant,"},
{"lineNum":"  277","line":"    bool& invertible"},
{"lineNum":"  278","line":"  )"},
{"lineNum":"  279","line":"  {"},
{"lineNum":"  280","line":"    using std::abs;"},
{"lineNum":"  281","line":"    determinant = matrix.determinant();"},
{"lineNum":"  282","line":"    invertible = abs(determinant) > absDeterminantThreshold;"},
{"lineNum":"  283","line":"    if(invertible && extract_data(matrix) != extract_data(inverse)) {"},
{"lineNum":"  284","line":"      compute_inverse<MatrixType, ResultType>::run(matrix, inverse);"},
{"lineNum":"  285","line":"    }"},
{"lineNum":"  286","line":"    else if(invertible) {"},
{"lineNum":"  287","line":"      MatrixType matrix_t = matrix;"},
{"lineNum":"  288","line":"      compute_inverse<MatrixType, ResultType>::run(matrix_t, inverse);"},
{"lineNum":"  289","line":"    }"},
{"lineNum":"  290","line":"  }"},
{"lineNum":"  291","line":"};"},
{"lineNum":"  292","line":""},
{"lineNum":"  293","line":"/*************************"},
{"lineNum":"  294","line":"*** MatrixBase methods ***"},
{"lineNum":"  295","line":"*************************/"},
{"lineNum":"  296","line":""},
{"lineNum":"  297","line":"} // end namespace internal"},
{"lineNum":"  298","line":""},
{"lineNum":"  299","line":"namespace internal {"},
{"lineNum":"  300","line":""},
{"lineNum":"  301","line":"// Specialization for \"dense = dense_xpr.inverse()\""},
{"lineNum":"  302","line":"template<typename DstXprType, typename XprType>"},
{"lineNum":"  303","line":"struct Assignment<DstXprType, Inverse<XprType>, internal::assign_op<typename DstXprType::Scalar,typename XprType::Scalar>, Dense2Dense>"},
{"lineNum":"  304","line":"{"},
{"lineNum":"  305","line":"  typedef Inverse<XprType> SrcXprType;"},
{"lineNum":"  306","line":"  EIGEN_DEVICE_FUNC"},
{"lineNum":"  307","line":"  static void run(DstXprType &dst, const SrcXprType &src, const internal::assign_op<typename DstXprType::Scalar,typename XprType::Scalar> &)"},
{"lineNum":"  308","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  309","line":"    Index dstRows = src.rows();"},
{"lineNum":"  310","line":"    Index dstCols = src.cols();"},
{"lineNum":"  311","line":"    if((dst.rows()!=dstRows) || (dst.cols()!=dstCols))","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  312","line":"      dst.resize(dstRows, dstCols);","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  313","line":""},
{"lineNum":"  314","line":"    const int Size = EIGEN_PLAIN_ENUM_MIN(XprType::ColsAtCompileTime,DstXprType::ColsAtCompileTime);"},
{"lineNum":"  315","line":"    EIGEN_ONLY_USED_FOR_DEBUG(Size);"},
{"lineNum":"  316","line":"    eigen_assert(( (Size<=1) || (Size>4) || (extract_data(src.nestedExpression())!=extract_data(dst)))"},
{"lineNum":"  317","line":"              && \"Aliasing problem detected in inverse(), you need to do inverse().eval() here.\");"},
{"lineNum":"  318","line":""},
{"lineNum":"  319","line":"    typedef typename internal::nested_eval<XprType,XprType::ColsAtCompileTime>::type  ActualXprType;"},
{"lineNum":"  320","line":"    typedef typename internal::remove_all<ActualXprType>::type                        ActualXprTypeCleanded;"},
{"lineNum":"  321","line":""},
{"lineNum":"  322","line":"    ActualXprType actual_xpr(src.nestedExpression());"},
{"lineNum":"  323","line":""},
{"lineNum":"  324","line":"    compute_inverse<ActualXprTypeCleanded, DstXprType>::run(actual_xpr, dst);","class":"lineNoCov","hits":"0","possible_hits":"6",},
{"lineNum":"  325","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  326","line":"};"},
{"lineNum":"  327","line":""},
{"lineNum":"  328","line":""},
{"lineNum":"  329","line":"} // end namespace internal"},
{"lineNum":"  330","line":""},
{"lineNum":"  331","line":"/** \\lu_module"},
{"lineNum":"  332","line":"  *"},
{"lineNum":"  333","line":"  * \\returns the matrix inverse of this matrix."},
{"lineNum":"  334","line":"  *"},
{"lineNum":"  335","line":"  * For small fixed sizes up to 4x4, this method uses cofactors."},
{"lineNum":"  336","line":"  * In the general case, this method uses class PartialPivLU."},
{"lineNum":"  337","line":"  *"},
{"lineNum":"  338","line":"  * \\note This matrix must be invertible, otherwise the result is undefined. If you need an"},
{"lineNum":"  339","line":"  * invertibility check, do the following:"},
{"lineNum":"  340","line":"  * \\li for fixed sizes up to 4x4, use computeInverseAndDetWithCheck()."},
{"lineNum":"  341","line":"  * \\li for the general case, use class FullPivLU."},
{"lineNum":"  342","line":"  *"},
{"lineNum":"  343","line":"  * Example: \\include MatrixBase_inverse.cpp"},
{"lineNum":"  344","line":"  * Output: \\verbinclude MatrixBase_inverse.out"},
{"lineNum":"  345","line":"  *"},
{"lineNum":"  346","line":"  * \\sa computeInverseAndDetWithCheck()"},
{"lineNum":"  347","line":"  */"},
{"lineNum":"  348","line":"template<typename Derived>"},
{"lineNum":"  349","line":"EIGEN_DEVICE_FUNC"},
{"lineNum":"  350","line":"inline const Inverse<Derived> MatrixBase<Derived>::inverse() const"},
{"lineNum":"  351","line":"{"},
{"lineNum":"  352","line":"  EIGEN_STATIC_ASSERT(!NumTraits<Scalar>::IsInteger,THIS_FUNCTION_IS_NOT_FOR_INTEGER_NUMERIC_TYPES)"},
{"lineNum":"  353","line":"  eigen_assert(rows() == cols());"},
{"lineNum":"  354","line":"  return Inverse<Derived>(derived());"},
{"lineNum":"  355","line":"}"},
{"lineNum":"  356","line":""},
{"lineNum":"  357","line":"/** \\lu_module"},
{"lineNum":"  358","line":"  *"},
{"lineNum":"  359","line":"  * Computation of matrix inverse and determinant, with invertibility check."},
{"lineNum":"  360","line":"  *"},
{"lineNum":"  361","line":"  * This is only for fixed-size square matrices of size up to 4x4."},
{"lineNum":"  362","line":"  *"},
{"lineNum":"  363","line":"  * Notice that it will trigger a copy of input matrix when trying to do the inverse in place."},
{"lineNum":"  364","line":"  *"},
{"lineNum":"  365","line":"  * \\param inverse Reference to the matrix in which to store the inverse."},
{"lineNum":"  366","line":"  * \\param determinant Reference to the variable in which to store the determinant."},
{"lineNum":"  367","line":"  * \\param invertible Reference to the bool variable in which to store whether the matrix is invertible."},
{"lineNum":"  368","line":"  * \\param absDeterminantThreshold Optional parameter controlling the invertibility check."},
{"lineNum":"  369","line":"  *                                The matrix will be declared invertible if the absolute value of its"},
{"lineNum":"  370","line":"  *                                determinant is greater than this threshold."},
{"lineNum":"  371","line":"  *"},
{"lineNum":"  372","line":"  * Example: \\include MatrixBase_computeInverseAndDetWithCheck.cpp"},
{"lineNum":"  373","line":"  * Output: \\verbinclude MatrixBase_computeInverseAndDetWithCheck.out"},
{"lineNum":"  374","line":"  *"},
{"lineNum":"  375","line":"  * \\sa inverse(), computeInverseWithCheck()"},
{"lineNum":"  376","line":"  */"},
{"lineNum":"  377","line":"template<typename Derived>"},
{"lineNum":"  378","line":"template<typename ResultType>"},
{"lineNum":"  379","line":"inline void MatrixBase<Derived>::computeInverseAndDetWithCheck("},
{"lineNum":"  380","line":"    ResultType& inverse,"},
{"lineNum":"  381","line":"    typename ResultType::Scalar& determinant,"},
{"lineNum":"  382","line":"    bool& invertible,"},
{"lineNum":"  383","line":"    const RealScalar& absDeterminantThreshold"},
{"lineNum":"  384","line":"  ) const"},
{"lineNum":"  385","line":"{"},
{"lineNum":"  386","line":"  // i\'d love to put some static assertions there, but SFINAE means that they have no effect..."},
{"lineNum":"  387","line":"  eigen_assert(rows() == cols());"},
{"lineNum":"  388","line":"  // for 2x2, it\'s worth giving a chance to avoid evaluating."},
{"lineNum":"  389","line":"  // for larger sizes, evaluating has negligible cost and limits code size."},
{"lineNum":"  390","line":"  typedef typename internal::conditional<"},
{"lineNum":"  391","line":"    RowsAtCompileTime == 2,"},
{"lineNum":"  392","line":"    typename internal::remove_all<typename internal::nested_eval<Derived, 2>::type>::type,"},
{"lineNum":"  393","line":"    PlainObject"},
{"lineNum":"  394","line":"  >::type MatrixType;"},
{"lineNum":"  395","line":"  internal::compute_inverse_and_det_with_check<MatrixType, ResultType>::run"},
{"lineNum":"  396","line":"    (derived(), absDeterminantThreshold, inverse, determinant, invertible);"},
{"lineNum":"  397","line":"}"},
{"lineNum":"  398","line":""},
{"lineNum":"  399","line":"/** \\lu_module"},
{"lineNum":"  400","line":"  *"},
{"lineNum":"  401","line":"  * Computation of matrix inverse, with invertibility check."},
{"lineNum":"  402","line":"  *"},
{"lineNum":"  403","line":"  * This is only for fixed-size square matrices of size up to 4x4."},
{"lineNum":"  404","line":"  *"},
{"lineNum":"  405","line":"  * Notice that it will trigger a copy of input matrix when trying to do the inverse in place."},
{"lineNum":"  406","line":"  *"},
{"lineNum":"  407","line":"  * \\param inverse Reference to the matrix in which to store the inverse."},
{"lineNum":"  408","line":"  * \\param invertible Reference to the bool variable in which to store whether the matrix is invertible."},
{"lineNum":"  409","line":"  * \\param absDeterminantThreshold Optional parameter controlling the invertibility check."},
{"lineNum":"  410","line":"  *                                The matrix will be declared invertible if the absolute value of its"},
{"lineNum":"  411","line":"  *                                determinant is greater than this threshold."},
{"lineNum":"  412","line":"  *"},
{"lineNum":"  413","line":"  * Example: \\include MatrixBase_computeInverseWithCheck.cpp"},
{"lineNum":"  414","line":"  * Output: \\verbinclude MatrixBase_computeInverseWithCheck.out"},
{"lineNum":"  415","line":"  *"},
{"lineNum":"  416","line":"  * \\sa inverse(), computeInverseAndDetWithCheck()"},
{"lineNum":"  417","line":"  */"},
{"lineNum":"  418","line":"template<typename Derived>"},
{"lineNum":"  419","line":"template<typename ResultType>"},
{"lineNum":"  420","line":"inline void MatrixBase<Derived>::computeInverseWithCheck("},
{"lineNum":"  421","line":"    ResultType& inverse,"},
{"lineNum":"  422","line":"    bool& invertible,"},
{"lineNum":"  423","line":"    const RealScalar& absDeterminantThreshold"},
{"lineNum":"  424","line":"  ) const"},
{"lineNum":"  425","line":"{"},
{"lineNum":"  426","line":"  Scalar determinant;"},
{"lineNum":"  427","line":"  // i\'d love to put some static assertions there, but SFINAE means that they have no effect..."},
{"lineNum":"  428","line":"  eigen_assert(rows() == cols());"},
{"lineNum":"  429","line":"  computeInverseAndDetWithCheck(inverse,determinant,invertible,absDeterminantThreshold);"},
{"lineNum":"  430","line":"}"},
{"lineNum":"  431","line":""},
{"lineNum":"  432","line":"} // end namespace Eigen"},
{"lineNum":"  433","line":""},
{"lineNum":"  434","line":"#endif // EIGEN_INVERSE_IMPL_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 8, "covered" : 0,};
var merged_data = [];
