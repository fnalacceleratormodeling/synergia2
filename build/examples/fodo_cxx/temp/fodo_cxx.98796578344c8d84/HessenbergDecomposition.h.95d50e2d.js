var data = {lines:[
{"lineNum":"    1","line":"// This file is part of Eigen, a lightweight C++ template library"},
{"lineNum":"    2","line":"// for linear algebra."},
{"lineNum":"    3","line":"//"},
{"lineNum":"    4","line":"// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>"},
{"lineNum":"    5","line":"// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>"},
{"lineNum":"    6","line":"//"},
{"lineNum":"    7","line":"// This Source Code Form is subject to the terms of the Mozilla"},
{"lineNum":"    8","line":"// Public License v. 2.0. If a copy of the MPL was not distributed"},
{"lineNum":"    9","line":"// with this file, You can obtain one at http://mozilla.org/MPL/2.0/."},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"#ifndef EIGEN_HESSENBERGDECOMPOSITION_H"},
{"lineNum":"   12","line":"#define EIGEN_HESSENBERGDECOMPOSITION_H"},
{"lineNum":"   13","line":""},
{"lineNum":"   14","line":"#include \"./InternalHeaderCheck.h\""},
{"lineNum":"   15","line":""},
{"lineNum":"   16","line":"namespace Eigen {"},
{"lineNum":"   17","line":""},
{"lineNum":"   18","line":"namespace internal {"},
{"lineNum":"   19","line":""},
{"lineNum":"   20","line":"template<typename MatrixType> struct HessenbergDecompositionMatrixHReturnType;"},
{"lineNum":"   21","line":"template<typename MatrixType>"},
{"lineNum":"   22","line":"struct traits<HessenbergDecompositionMatrixHReturnType<MatrixType> >"},
{"lineNum":"   23","line":"{"},
{"lineNum":"   24","line":"  typedef MatrixType ReturnType;"},
{"lineNum":"   25","line":"};"},
{"lineNum":"   26","line":""},
{"lineNum":"   27","line":"}"},
{"lineNum":"   28","line":""},
{"lineNum":"   29","line":"/** \\eigenvalues_module \\ingroup Eigenvalues_Module"},
{"lineNum":"   30","line":"  *"},
{"lineNum":"   31","line":"  *"},
{"lineNum":"   32","line":"  * \\class HessenbergDecomposition"},
{"lineNum":"   33","line":"  *"},
{"lineNum":"   34","line":"  * \\brief Reduces a square matrix to Hessenberg form by an orthogonal similarity transformation"},
{"lineNum":"   35","line":"  *"},
{"lineNum":"   36","line":"  * \\tparam MatrixType_ the type of the matrix of which we are computing the Hessenberg decomposition"},
{"lineNum":"   37","line":"  *"},
{"lineNum":"   38","line":"  * This class performs an Hessenberg decomposition of a matrix \\f$ A \\f$. In"},
{"lineNum":"   39","line":"  * the real case, the Hessenberg decomposition consists of an orthogonal"},
{"lineNum":"   40","line":"  * matrix \\f$ Q \\f$ and a Hessenberg matrix \\f$ H \\f$ such that \\f$ A = Q H"},
{"lineNum":"   41","line":"  * Q^T \\f$. An orthogonal matrix is a matrix whose inverse equals its"},
{"lineNum":"   42","line":"  * transpose (\\f$ Q^{-1} = Q^T \\f$). A Hessenberg matrix has zeros below the"},
{"lineNum":"   43","line":"  * subdiagonal, so it is almost upper triangular. The Hessenberg decomposition"},
{"lineNum":"   44","line":"  * of a complex matrix is \\f$ A = Q H Q^* \\f$ with \\f$ Q \\f$ unitary (that is,"},
{"lineNum":"   45","line":"  * \\f$ Q^{-1} = Q^* \\f$)."},
{"lineNum":"   46","line":"  *"},
{"lineNum":"   47","line":"  * Call the function compute() to compute the Hessenberg decomposition of a"},
{"lineNum":"   48","line":"  * given matrix. Alternatively, you can use the"},
{"lineNum":"   49","line":"  * HessenbergDecomposition(const MatrixType&) constructor which computes the"},
{"lineNum":"   50","line":"  * Hessenberg decomposition at construction time. Once the decomposition is"},
{"lineNum":"   51","line":"  * computed, you can use the matrixH() and matrixQ() functions to construct"},
{"lineNum":"   52","line":"  * the matrices H and Q in the decomposition."},
{"lineNum":"   53","line":"  *"},
{"lineNum":"   54","line":"  * The documentation for matrixH() contains an example of the typical use of"},
{"lineNum":"   55","line":"  * this class."},
{"lineNum":"   56","line":"  *"},
{"lineNum":"   57","line":"  * \\sa class ComplexSchur, class Tridiagonalization, \\ref QR_Module \"QR Module\""},
{"lineNum":"   58","line":"  */"},
{"lineNum":"   59","line":"template<typename MatrixType_> class HessenbergDecomposition"},
{"lineNum":"   60","line":"{"},
{"lineNum":"   61","line":"  public:"},
{"lineNum":"   62","line":""},
{"lineNum":"   63","line":"    /** \\brief Synonym for the template parameter \\p MatrixType_. */"},
{"lineNum":"   64","line":"    typedef MatrixType_ MatrixType;"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"    enum {"},
{"lineNum":"   67","line":"      Size = MatrixType::RowsAtCompileTime,"},
{"lineNum":"   68","line":"      SizeMinusOne = Size == Dynamic ? Dynamic : Size - 1,"},
{"lineNum":"   69","line":"      Options = MatrixType::Options,"},
{"lineNum":"   70","line":"      MaxSize = MatrixType::MaxRowsAtCompileTime,"},
{"lineNum":"   71","line":"      MaxSizeMinusOne = MaxSize == Dynamic ? Dynamic : MaxSize - 1"},
{"lineNum":"   72","line":"    };"},
{"lineNum":"   73","line":""},
{"lineNum":"   74","line":"    /** \\brief Scalar type for matrices of type #MatrixType. */"},
{"lineNum":"   75","line":"    typedef typename MatrixType::Scalar Scalar;"},
{"lineNum":"   76","line":"    typedef Eigen::Index Index; ///< \\deprecated since Eigen 3.3"},
{"lineNum":"   77","line":""},
{"lineNum":"   78","line":"    /** \\brief Type for vector of Householder coefficients."},
{"lineNum":"   79","line":"      *"},
{"lineNum":"   80","line":"      * This is column vector with entries of type #Scalar. The length of the"},
{"lineNum":"   81","line":"      * vector is one less than the size of #MatrixType, if it is a fixed-side"},
{"lineNum":"   82","line":"      * type."},
{"lineNum":"   83","line":"      */"},
{"lineNum":"   84","line":"    typedef Matrix<Scalar, SizeMinusOne, 1, Options & ~RowMajor, MaxSizeMinusOne, 1> CoeffVectorType;"},
{"lineNum":"   85","line":""},
{"lineNum":"   86","line":"    /** \\brief Return type of matrixQ() */"},
{"lineNum":"   87","line":"    typedef HouseholderSequence<MatrixType,typename internal::remove_all<typename CoeffVectorType::ConjugateReturnType>::type> HouseholderSequenceType;"},
{"lineNum":"   88","line":""},
{"lineNum":"   89","line":"    typedef internal::HessenbergDecompositionMatrixHReturnType<MatrixType> MatrixHReturnType;"},
{"lineNum":"   90","line":""},
{"lineNum":"   91","line":"    /** \\brief Default constructor; the decomposition will be computed later."},
{"lineNum":"   92","line":"      *"},
{"lineNum":"   93","line":"      * \\param [in] size  The size of the matrix whose Hessenberg decomposition will be computed."},
{"lineNum":"   94","line":"      *"},
{"lineNum":"   95","line":"      * The default constructor is useful in cases in which the user intends to"},
{"lineNum":"   96","line":"      * perform decompositions via compute().  The \\p size parameter is only"},
{"lineNum":"   97","line":"      * used as a hint. It is not an error to give a wrong \\p size, but it may"},
{"lineNum":"   98","line":"      * impair performance."},
{"lineNum":"   99","line":"      *"},
{"lineNum":"  100","line":"      * \\sa compute() for an example."},
{"lineNum":"  101","line":"      */"},
{"lineNum":"  102","line":"    explicit HessenbergDecomposition(Index size = Size==Dynamic ? 2 : Size)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  103","line":"      : m_matrix(size,size),"},
{"lineNum":"  104","line":"        m_temp(size),","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  105","line":"        m_isInitialized(false)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  106","line":"    {","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  107","line":"      if(size>1)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  108","line":"        m_hCoeffs.resize(size-1);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  109","line":"    }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  110","line":""},
{"lineNum":"  111","line":"    /** \\brief Constructor; computes Hessenberg decomposition of given matrix."},
{"lineNum":"  112","line":"      *"},
{"lineNum":"  113","line":"      * \\param[in]  matrix  Square matrix whose Hessenberg decomposition is to be computed."},
{"lineNum":"  114","line":"      *"},
{"lineNum":"  115","line":"      * This constructor calls compute() to compute the Hessenberg"},
{"lineNum":"  116","line":"      * decomposition."},
{"lineNum":"  117","line":"      *"},
{"lineNum":"  118","line":"      * \\sa matrixH() for an example."},
{"lineNum":"  119","line":"      */"},
{"lineNum":"  120","line":"    template<typename InputType>"},
{"lineNum":"  121","line":"    explicit HessenbergDecomposition(const EigenBase<InputType>& matrix)"},
{"lineNum":"  122","line":"      : m_matrix(matrix.derived()),"},
{"lineNum":"  123","line":"        m_temp(matrix.rows()),"},
{"lineNum":"  124","line":"        m_isInitialized(false)"},
{"lineNum":"  125","line":"    {"},
{"lineNum":"  126","line":"      if(matrix.rows()<2)"},
{"lineNum":"  127","line":"      {"},
{"lineNum":"  128","line":"        m_isInitialized = true;"},
{"lineNum":"  129","line":"        return;"},
{"lineNum":"  130","line":"      }"},
{"lineNum":"  131","line":"      m_hCoeffs.resize(matrix.rows()-1,1);"},
{"lineNum":"  132","line":"      _compute(m_matrix, m_hCoeffs, m_temp);"},
{"lineNum":"  133","line":"      m_isInitialized = true;"},
{"lineNum":"  134","line":"    }"},
{"lineNum":"  135","line":""},
{"lineNum":"  136","line":"    /** \\brief Computes Hessenberg decomposition of given matrix."},
{"lineNum":"  137","line":"      *"},
{"lineNum":"  138","line":"      * \\param[in]  matrix  Square matrix whose Hessenberg decomposition is to be computed."},
{"lineNum":"  139","line":"      * \\returns    Reference to \\c *this"},
{"lineNum":"  140","line":"      *"},
{"lineNum":"  141","line":"      * The Hessenberg decomposition is computed by bringing the columns of the"},
{"lineNum":"  142","line":"      * matrix successively in the required form using Householder reflections"},
{"lineNum":"  143","line":"      * (see, e.g., Algorithm 7.4.2 in Golub \\& Van Loan, <i>%Matrix"},
{"lineNum":"  144","line":"      * Computations</i>). The cost is \\f$ 10n^3/3 \\f$ flops, where \\f$ n \\f$"},
{"lineNum":"  145","line":"      * denotes the size of the given matrix."},
{"lineNum":"  146","line":"      *"},
{"lineNum":"  147","line":"      * This method reuses of the allocated data in the HessenbergDecomposition"},
{"lineNum":"  148","line":"      * object."},
{"lineNum":"  149","line":"      *"},
{"lineNum":"  150","line":"      * Example: \\include HessenbergDecomposition_compute.cpp"},
{"lineNum":"  151","line":"      * Output: \\verbinclude HessenbergDecomposition_compute.out"},
{"lineNum":"  152","line":"      */"},
{"lineNum":"  153","line":"    template<typename InputType>"},
{"lineNum":"  154","line":"    HessenbergDecomposition& compute(const EigenBase<InputType>& matrix)"},
{"lineNum":"  155","line":"    {"},
{"lineNum":"  156","line":"      m_matrix = matrix.derived();"},
{"lineNum":"  157","line":"      if(matrix.rows()<2)"},
{"lineNum":"  158","line":"      {"},
{"lineNum":"  159","line":"        m_isInitialized = true;"},
{"lineNum":"  160","line":"        return *this;"},
{"lineNum":"  161","line":"      }"},
{"lineNum":"  162","line":"      m_hCoeffs.resize(matrix.rows()-1,1);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  163","line":"      _compute(m_matrix, m_hCoeffs, m_temp);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  164","line":"      m_isInitialized = true;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  165","line":"      return *this;"},
{"lineNum":"  166","line":"    }"},
{"lineNum":"  167","line":""},
{"lineNum":"  168","line":"    /** \\brief Returns the Householder coefficients."},
{"lineNum":"  169","line":"      *"},
{"lineNum":"  170","line":"      * \\returns a const reference to the vector of Householder coefficients"},
{"lineNum":"  171","line":"      *"},
{"lineNum":"  172","line":"      * \\pre Either the constructor HessenbergDecomposition(const MatrixType&)"},
{"lineNum":"  173","line":"      * or the member function compute(const MatrixType&) has been called"},
{"lineNum":"  174","line":"      * before to compute the Hessenberg decomposition of a matrix."},
{"lineNum":"  175","line":"      *"},
{"lineNum":"  176","line":"      * The Householder coefficients allow the reconstruction of the matrix"},
{"lineNum":"  177","line":"      * \\f$ Q \\f$ in the Hessenberg decomposition from the packed data."},
{"lineNum":"  178","line":"      *"},
{"lineNum":"  179","line":"      * \\sa packedMatrix(), \\ref Householder_Module \"Householder module\""},
{"lineNum":"  180","line":"      */"},
{"lineNum":"  181","line":"    const CoeffVectorType& householderCoefficients() const"},
{"lineNum":"  182","line":"    {"},
{"lineNum":"  183","line":"      eigen_assert(m_isInitialized && \"HessenbergDecomposition is not initialized.\");"},
{"lineNum":"  184","line":"      return m_hCoeffs;"},
{"lineNum":"  185","line":"    }"},
{"lineNum":"  186","line":""},
{"lineNum":"  187","line":"    /** \\brief Returns the internal representation of the decomposition"},
{"lineNum":"  188","line":"      *"},
{"lineNum":"  189","line":"      *\t\\returns a const reference to a matrix with the internal representation"},
{"lineNum":"  190","line":"      *\t         of the decomposition."},
{"lineNum":"  191","line":"      *"},
{"lineNum":"  192","line":"      * \\pre Either the constructor HessenbergDecomposition(const MatrixType&)"},
{"lineNum":"  193","line":"      * or the member function compute(const MatrixType&) has been called"},
{"lineNum":"  194","line":"      * before to compute the Hessenberg decomposition of a matrix."},
{"lineNum":"  195","line":"      *"},
{"lineNum":"  196","line":"      * The returned matrix contains the following information:"},
{"lineNum":"  197","line":"      *  - the upper part and lower sub-diagonal represent the Hessenberg matrix H"},
{"lineNum":"  198","line":"      *  - the rest of the lower part contains the Householder vectors that, combined with"},
{"lineNum":"  199","line":"      *    Householder coefficients returned by householderCoefficients(),"},
{"lineNum":"  200","line":"      *    allows to reconstruct the matrix Q as"},
{"lineNum":"  201","line":"      *       \\f$ Q = H_{N-1} \\ldots H_1 H_0 \\f$."},
{"lineNum":"  202","line":"      *    Here, the matrices \\f$ H_i \\f$ are the Householder transformations"},
{"lineNum":"  203","line":"      *       \\f$ H_i = (I - h_i v_i v_i^T) \\f$"},
{"lineNum":"  204","line":"      *    where \\f$ h_i \\f$ is the \\f$ i \\f$th Householder coefficient and"},
{"lineNum":"  205","line":"      *    \\f$ v_i \\f$ is the Householder vector defined by"},
{"lineNum":"  206","line":"      *       \\f$ v_i = [ 0, \\ldots, 0, 1, M(i+2,i), \\ldots, M(N-1,i) ]^T \\f$"},
{"lineNum":"  207","line":"      *    with M the matrix returned by this function."},
{"lineNum":"  208","line":"      *"},
{"lineNum":"  209","line":"      * See LAPACK for further details on this packed storage."},
{"lineNum":"  210","line":"      *"},
{"lineNum":"  211","line":"      * Example: \\include HessenbergDecomposition_packedMatrix.cpp"},
{"lineNum":"  212","line":"      * Output: \\verbinclude HessenbergDecomposition_packedMatrix.out"},
{"lineNum":"  213","line":"      *"},
{"lineNum":"  214","line":"      * \\sa householderCoefficients()"},
{"lineNum":"  215","line":"      */"},
{"lineNum":"  216","line":"    const MatrixType& packedMatrix() const"},
{"lineNum":"  217","line":"    {"},
{"lineNum":"  218","line":"      eigen_assert(m_isInitialized && \"HessenbergDecomposition is not initialized.\");"},
{"lineNum":"  219","line":"      return m_matrix;"},
{"lineNum":"  220","line":"    }"},
{"lineNum":"  221","line":""},
{"lineNum":"  222","line":"    /** \\brief Reconstructs the orthogonal matrix Q in the decomposition"},
{"lineNum":"  223","line":"      *"},
{"lineNum":"  224","line":"      * \\returns object representing the matrix Q"},
{"lineNum":"  225","line":"      *"},
{"lineNum":"  226","line":"      * \\pre Either the constructor HessenbergDecomposition(const MatrixType&)"},
{"lineNum":"  227","line":"      * or the member function compute(const MatrixType&) has been called"},
{"lineNum":"  228","line":"      * before to compute the Hessenberg decomposition of a matrix."},
{"lineNum":"  229","line":"      *"},
{"lineNum":"  230","line":"      * This function returns a light-weight object of template class"},
{"lineNum":"  231","line":"      * HouseholderSequence. You can either apply it directly to a matrix or"},
{"lineNum":"  232","line":"      * you can convert it to a matrix of type #MatrixType."},
{"lineNum":"  233","line":"      *"},
{"lineNum":"  234","line":"      * \\sa matrixH() for an example, class HouseholderSequence"},
{"lineNum":"  235","line":"      */"},
{"lineNum":"  236","line":"    HouseholderSequenceType matrixQ() const"},
{"lineNum":"  237","line":"    {"},
{"lineNum":"  238","line":"      eigen_assert(m_isInitialized && \"HessenbergDecomposition is not initialized.\");"},
{"lineNum":"  239","line":"      return HouseholderSequenceType(m_matrix, m_hCoeffs.conjugate())"},
{"lineNum":"  240","line":"             .setLength(m_matrix.rows() - 1)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  241","line":"             .setShift(1);"},
{"lineNum":"  242","line":"    }"},
{"lineNum":"  243","line":""},
{"lineNum":"  244","line":"    /** \\brief Constructs the Hessenberg matrix H in the decomposition"},
{"lineNum":"  245","line":"      *"},
{"lineNum":"  246","line":"      * \\returns expression object representing the matrix H"},
{"lineNum":"  247","line":"      *"},
{"lineNum":"  248","line":"      * \\pre Either the constructor HessenbergDecomposition(const MatrixType&)"},
{"lineNum":"  249","line":"      * or the member function compute(const MatrixType&) has been called"},
{"lineNum":"  250","line":"      * before to compute the Hessenberg decomposition of a matrix."},
{"lineNum":"  251","line":"      *"},
{"lineNum":"  252","line":"      * The object returned by this function constructs the Hessenberg matrix H"},
{"lineNum":"  253","line":"      * when it is assigned to a matrix or otherwise evaluated. The matrix H is"},
{"lineNum":"  254","line":"      * constructed from the packed matrix as returned by packedMatrix(): The"},
{"lineNum":"  255","line":"      * upper part (including the subdiagonal) of the packed matrix contains"},
{"lineNum":"  256","line":"      * the matrix H. It may sometimes be better to directly use the packed"},
{"lineNum":"  257","line":"      * matrix instead of constructing the matrix H."},
{"lineNum":"  258","line":"      *"},
{"lineNum":"  259","line":"      * Example: \\include HessenbergDecomposition_matrixH.cpp"},
{"lineNum":"  260","line":"      * Output: \\verbinclude HessenbergDecomposition_matrixH.out"},
{"lineNum":"  261","line":"      *"},
{"lineNum":"  262","line":"      * \\sa matrixQ(), packedMatrix()"},
{"lineNum":"  263","line":"      */"},
{"lineNum":"  264","line":"    MatrixHReturnType matrixH() const"},
{"lineNum":"  265","line":"    {"},
{"lineNum":"  266","line":"      eigen_assert(m_isInitialized && \"HessenbergDecomposition is not initialized.\");"},
{"lineNum":"  267","line":"      return MatrixHReturnType(*this);"},
{"lineNum":"  268","line":"    }"},
{"lineNum":"  269","line":""},
{"lineNum":"  270","line":"  private:"},
{"lineNum":"  271","line":""},
{"lineNum":"  272","line":"    typedef Matrix<Scalar, 1, Size, int(Options) | int(RowMajor), 1, MaxSize> VectorType;"},
{"lineNum":"  273","line":"    typedef typename NumTraits<Scalar>::Real RealScalar;"},
{"lineNum":"  274","line":"    static void _compute(MatrixType& matA, CoeffVectorType& hCoeffs, VectorType& temp);"},
{"lineNum":"  275","line":""},
{"lineNum":"  276","line":"  protected:"},
{"lineNum":"  277","line":"    MatrixType m_matrix;"},
{"lineNum":"  278","line":"    CoeffVectorType m_hCoeffs;"},
{"lineNum":"  279","line":"    VectorType m_temp;"},
{"lineNum":"  280","line":"    bool m_isInitialized;"},
{"lineNum":"  281","line":"};"},
{"lineNum":"  282","line":""},
{"lineNum":"  283","line":"/** \\internal"},
{"lineNum":"  284","line":"  * Performs a tridiagonal decomposition of \\a matA in place."},
{"lineNum":"  285","line":"  *"},
{"lineNum":"  286","line":"  * \\param matA the input selfadjoint matrix"},
{"lineNum":"  287","line":"  * \\param hCoeffs returned Householder coefficients"},
{"lineNum":"  288","line":"  *"},
{"lineNum":"  289","line":"  * The result is written in the lower triangular part of \\a matA."},
{"lineNum":"  290","line":"  *"},
{"lineNum":"  291","line":"  * Implemented from Golub\'s \"%Matrix Computations\", algorithm 8.3.1."},
{"lineNum":"  292","line":"  *"},
{"lineNum":"  293","line":"  * \\sa packedMatrix()"},
{"lineNum":"  294","line":"  */"},
{"lineNum":"  295","line":"template<typename MatrixType>"},
{"lineNum":"  296","line":"void HessenbergDecomposition<MatrixType>::_compute(MatrixType& matA, CoeffVectorType& hCoeffs, VectorType& temp)"},
{"lineNum":"  297","line":"{","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  298","line":"  eigen_assert(matA.rows()==matA.cols());"},
{"lineNum":"  299","line":"  Index n = matA.rows();"},
{"lineNum":"  300","line":"  temp.resize(n);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  301","line":"  for (Index i = 0; i<n-1; ++i)","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  302","line":"  {"},
{"lineNum":"  303","line":"    // let\'s consider the vector v = i-th column starting at position i+1"},
{"lineNum":"  304","line":"    Index remainingSize = n-i-1;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  305","line":"    RealScalar beta;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  306","line":"    Scalar h;"},
{"lineNum":"  307","line":"    matA.col(i).tail(remainingSize).makeHouseholderInPlace(h, beta);"},
{"lineNum":"  308","line":"    matA.col(i).coeffRef(i+1) = beta;","class":"lineNoCov","hits":"0","possible_hits":"5",},
{"lineNum":"  309","line":"    hCoeffs.coeffRef(i) = h;","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  310","line":""},
{"lineNum":"  311","line":"    // Apply similarity transformation to remaining columns,"},
{"lineNum":"  312","line":"    // i.e., compute A = H A H\'"},
{"lineNum":"  313","line":""},
{"lineNum":"  314","line":"    // A = H A"},
{"lineNum":"  315","line":"    matA.bottomRightCorner(remainingSize, remainingSize)"},
{"lineNum":"  316","line":"        .applyHouseholderOnTheLeft(matA.col(i).tail(remainingSize-1), h, &temp.coeffRef(0));","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  317","line":""},
{"lineNum":"  318","line":"    // A = A H\'"},
{"lineNum":"  319","line":"    matA.rightCols(remainingSize)"},
{"lineNum":"  320","line":"        .applyHouseholderOnTheRight(matA.col(i).tail(remainingSize-1), numext::conj(h), &temp.coeffRef(0));","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  321","line":"  }"},
{"lineNum":"  322","line":"}","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  323","line":""},
{"lineNum":"  324","line":"namespace internal {"},
{"lineNum":"  325","line":""},
{"lineNum":"  326","line":"/** \\eigenvalues_module \\ingroup Eigenvalues_Module"},
{"lineNum":"  327","line":"  *"},
{"lineNum":"  328","line":"  *"},
{"lineNum":"  329","line":"  * \\brief Expression type for return value of HessenbergDecomposition::matrixH()"},
{"lineNum":"  330","line":"  *"},
{"lineNum":"  331","line":"  * \\tparam MatrixType type of matrix in the Hessenberg decomposition"},
{"lineNum":"  332","line":"  *"},
{"lineNum":"  333","line":"  * Objects of this type represent the Hessenberg matrix in the Hessenberg"},
{"lineNum":"  334","line":"  * decomposition of some matrix. The object holds a reference to the"},
{"lineNum":"  335","line":"  * HessenbergDecomposition class until the it is assigned or evaluated for"},
{"lineNum":"  336","line":"  * some other reason (the reference should remain valid during the life time"},
{"lineNum":"  337","line":"  * of this object). This class is the return type of"},
{"lineNum":"  338","line":"  * HessenbergDecomposition::matrixH(); there is probably no other use for this"},
{"lineNum":"  339","line":"  * class."},
{"lineNum":"  340","line":"  */"},
{"lineNum":"  341","line":"template<typename MatrixType> struct HessenbergDecompositionMatrixHReturnType"},
{"lineNum":"  342","line":": public ReturnByValue<HessenbergDecompositionMatrixHReturnType<MatrixType> >"},
{"lineNum":"  343","line":"{"},
{"lineNum":"  344","line":"  public:"},
{"lineNum":"  345","line":"    /** \\brief Constructor."},
{"lineNum":"  346","line":"      *"},
{"lineNum":"  347","line":"      * \\param[in] hess  Hessenberg decomposition"},
{"lineNum":"  348","line":"      */"},
{"lineNum":"  349","line":"    HessenbergDecompositionMatrixHReturnType(const HessenbergDecomposition<MatrixType>& hess) : m_hess(hess) { }"},
{"lineNum":"  350","line":""},
{"lineNum":"  351","line":"    /** \\brief Hessenberg matrix in decomposition."},
{"lineNum":"  352","line":"      *"},
{"lineNum":"  353","line":"      * \\param[out] result  Hessenberg matrix in decomposition \\p hess which"},
{"lineNum":"  354","line":"      *                     was passed to the constructor"},
{"lineNum":"  355","line":"      */"},
{"lineNum":"  356","line":"    template <typename ResultType>"},
{"lineNum":"  357","line":"    inline void evalTo(ResultType& result) const"},
{"lineNum":"  358","line":"    {"},
{"lineNum":"  359","line":"      result = m_hess.packedMatrix();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  360","line":"      Index n = result.rows();"},
{"lineNum":"  361","line":"      if (n>2)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  362","line":"        result.bottomLeftCorner(n-2, n-2).template triangularView<Lower>().setZero();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  363","line":"    }"},
{"lineNum":"  364","line":""},
{"lineNum":"  365","line":"    Index rows() const { return m_hess.packedMatrix().rows(); }","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  366","line":"    Index cols() const { return m_hess.packedMatrix().cols(); }"},
{"lineNum":"  367","line":""},
{"lineNum":"  368","line":"  protected:"},
{"lineNum":"  369","line":"    const HessenbergDecomposition<MatrixType>& m_hess;"},
{"lineNum":"  370","line":"};"},
{"lineNum":"  371","line":""},
{"lineNum":"  372","line":"} // end namespace internal"},
{"lineNum":"  373","line":""},
{"lineNum":"  374","line":"} // end namespace Eigen"},
{"lineNum":"  375","line":""},
{"lineNum":"  376","line":"#endif // EIGEN_HESSENBERGDECOMPOSITION_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 25, "covered" : 0,};
var merged_data = [];
