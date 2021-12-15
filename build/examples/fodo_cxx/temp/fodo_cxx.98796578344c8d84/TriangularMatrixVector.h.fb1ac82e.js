var data = {lines:[
{"lineNum":"    1","line":"// This file is part of Eigen, a lightweight C++ template library"},
{"lineNum":"    2","line":"// for linear algebra."},
{"lineNum":"    3","line":"//"},
{"lineNum":"    4","line":"// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>"},
{"lineNum":"    5","line":"//"},
{"lineNum":"    6","line":"// This Source Code Form is subject to the terms of the Mozilla"},
{"lineNum":"    7","line":"// Public License v. 2.0. If a copy of the MPL was not distributed"},
{"lineNum":"    8","line":"// with this file, You can obtain one at http://mozilla.org/MPL/2.0/."},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"#ifndef EIGEN_TRIANGULARMATRIXVECTOR_H"},
{"lineNum":"   11","line":"#define EIGEN_TRIANGULARMATRIXVECTOR_H"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#include \"../InternalHeaderCheck.h\""},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"namespace Eigen {"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"namespace internal {"},
{"lineNum":"   18","line":""},
{"lineNum":"   19","line":"template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs, int StorageOrder, int Version=Specialized>"},
{"lineNum":"   20","line":"struct triangular_matrix_vector_product;"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs, int Version>"},
{"lineNum":"   23","line":"struct triangular_matrix_vector_product<Index,Mode,LhsScalar,ConjLhs,RhsScalar,ConjRhs,ColMajor,Version>"},
{"lineNum":"   24","line":"{"},
{"lineNum":"   25","line":"  typedef typename ScalarBinaryOpTraits<LhsScalar, RhsScalar>::ReturnType ResScalar;"},
{"lineNum":"   26","line":"  enum {"},
{"lineNum":"   27","line":"    IsLower = ((Mode&Lower)==Lower),"},
{"lineNum":"   28","line":"    HasUnitDiag = (Mode & UnitDiag)==UnitDiag,"},
{"lineNum":"   29","line":"    HasZeroDiag = (Mode & ZeroDiag)==ZeroDiag"},
{"lineNum":"   30","line":"  };"},
{"lineNum":"   31","line":"  static EIGEN_DONT_INLINE  void run(Index _rows, Index _cols, const LhsScalar* _lhs, Index lhsStride,"},
{"lineNum":"   32","line":"                                     const RhsScalar* _rhs, Index rhsIncr, ResScalar* _res, Index resIncr, const RhsScalar& alpha);"},
{"lineNum":"   33","line":"};"},
{"lineNum":"   34","line":""},
{"lineNum":"   35","line":"template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs, int Version>"},
{"lineNum":"   36","line":"EIGEN_DONT_INLINE void triangular_matrix_vector_product<Index,Mode,LhsScalar,ConjLhs,RhsScalar,ConjRhs,ColMajor,Version>"},
{"lineNum":"   37","line":"  ::run(Index _rows, Index _cols, const LhsScalar* _lhs, Index lhsStride,"},
{"lineNum":"   38","line":"        const RhsScalar* _rhs, Index rhsIncr, ResScalar* _res, Index resIncr, const RhsScalar& alpha)"},
{"lineNum":"   39","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   40","line":"    static const Index PanelWidth = EIGEN_TUNE_TRIANGULAR_PANEL_WIDTH;"},
{"lineNum":"   41","line":"    Index size = (std::min)(_rows,_cols);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   42","line":"    Index rows = IsLower ? _rows : (std::min)(_rows,_cols);"},
{"lineNum":"   43","line":"    Index cols = IsLower ? (std::min)(_rows,_cols) : _cols;"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"    typedef Map<const Matrix<LhsScalar,Dynamic,Dynamic,ColMajor>, 0, OuterStride<> > LhsMap;"},
{"lineNum":"   46","line":"    const LhsMap lhs(_lhs,rows,cols,OuterStride<>(lhsStride));"},
{"lineNum":"   47","line":"    typename conj_expr_if<ConjLhs,LhsMap>::type cjLhs(lhs);"},
{"lineNum":"   48","line":""},
{"lineNum":"   49","line":"    typedef Map<const Matrix<RhsScalar,Dynamic,1>, 0, InnerStride<> > RhsMap;"},
{"lineNum":"   50","line":"    const RhsMap rhs(_rhs,cols,InnerStride<>(rhsIncr));"},
{"lineNum":"   51","line":"    typename conj_expr_if<ConjRhs,RhsMap>::type cjRhs(rhs);"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"    typedef Map<Matrix<ResScalar,Dynamic,1> > ResMap;"},
{"lineNum":"   54","line":"    ResMap res(_res,rows);"},
{"lineNum":"   55","line":""},
{"lineNum":"   56","line":"    typedef const_blas_data_mapper<LhsScalar,Index,ColMajor> LhsMapper;"},
{"lineNum":"   57","line":"    typedef const_blas_data_mapper<RhsScalar,Index,RowMajor> RhsMapper;"},
{"lineNum":"   58","line":""},
{"lineNum":"   59","line":"    for (Index pi=0; pi<size; pi+=PanelWidth)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   60","line":"    {"},
{"lineNum":"   61","line":"      Index actualPanelWidth = (std::min)(PanelWidth, size-pi);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   62","line":"      for (Index k=0; k<actualPanelWidth; ++k)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   63","line":"      {"},
{"lineNum":"   64","line":"        Index i = pi + k;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   65","line":"        Index s = IsLower ? ((HasUnitDiag||HasZeroDiag) ? i+1 : i ) : pi;"},
{"lineNum":"   66","line":"        Index r = IsLower ? actualPanelWidth-k : k+1;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   67","line":"        if ((!(HasUnitDiag||HasZeroDiag)) || (--r)>0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   68","line":"          res.segment(s,r) += (alpha * cjRhs.coeff(i)) * cjLhs.col(i).segment(s,r);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   69","line":"        if (HasUnitDiag)"},
{"lineNum":"   70","line":"          res.coeffRef(i) += alpha * cjRhs.coeff(i);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   71","line":"      }"},
{"lineNum":"   72","line":"      Index r = IsLower ? rows - pi - actualPanelWidth : pi;"},
{"lineNum":"   73","line":"      if (r>0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   74","line":"      {"},
{"lineNum":"   75","line":"        Index s = IsLower ? pi+actualPanelWidth : 0;"},
{"lineNum":"   76","line":"        general_matrix_vector_product<Index,LhsScalar,LhsMapper,ColMajor,ConjLhs,RhsScalar,RhsMapper,ConjRhs,BuiltIn>::run(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   77","line":"            r, actualPanelWidth,"},
{"lineNum":"   78","line":"            LhsMapper(&lhs.coeffRef(s,pi), lhsStride),"},
{"lineNum":"   79","line":"            RhsMapper(&rhs.coeffRef(pi), rhsIncr),"},
{"lineNum":"   80","line":"            &res.coeffRef(s), resIncr, alpha);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   81","line":"      }"},
{"lineNum":"   82","line":"    }"},
{"lineNum":"   83","line":"    if((!IsLower) && cols>size)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   84","line":"    {"},
{"lineNum":"   85","line":"      general_matrix_vector_product<Index,LhsScalar,LhsMapper,ColMajor,ConjLhs,RhsScalar,RhsMapper,ConjRhs>::run(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   86","line":"          rows, cols-size,","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   87","line":"          LhsMapper(&lhs.coeffRef(0,size), lhsStride),"},
{"lineNum":"   88","line":"          RhsMapper(&rhs.coeffRef(size), rhsIncr),"},
{"lineNum":"   89","line":"          _res, resIncr, alpha);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   90","line":"    }"},
{"lineNum":"   91","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"   92","line":""},
{"lineNum":"   93","line":"template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs,int Version>"},
{"lineNum":"   94","line":"struct triangular_matrix_vector_product<Index,Mode,LhsScalar,ConjLhs,RhsScalar,ConjRhs,RowMajor,Version>"},
{"lineNum":"   95","line":"{"},
{"lineNum":"   96","line":"  typedef typename ScalarBinaryOpTraits<LhsScalar, RhsScalar>::ReturnType ResScalar;"},
{"lineNum":"   97","line":"  enum {"},
{"lineNum":"   98","line":"    IsLower = ((Mode&Lower)==Lower),"},
{"lineNum":"   99","line":"    HasUnitDiag = (Mode & UnitDiag)==UnitDiag,"},
{"lineNum":"  100","line":"    HasZeroDiag = (Mode & ZeroDiag)==ZeroDiag"},
{"lineNum":"  101","line":"  };"},
{"lineNum":"  102","line":"  static EIGEN_DONT_INLINE void run(Index _rows, Index _cols, const LhsScalar* _lhs, Index lhsStride,"},
{"lineNum":"  103","line":"                                    const RhsScalar* _rhs, Index rhsIncr, ResScalar* _res, Index resIncr, const ResScalar& alpha);"},
{"lineNum":"  104","line":"};"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"template<typename Index, int Mode, typename LhsScalar, bool ConjLhs, typename RhsScalar, bool ConjRhs,int Version>"},
{"lineNum":"  107","line":"EIGEN_DONT_INLINE void triangular_matrix_vector_product<Index,Mode,LhsScalar,ConjLhs,RhsScalar,ConjRhs,RowMajor,Version>"},
{"lineNum":"  108","line":"  ::run(Index _rows, Index _cols, const LhsScalar* _lhs, Index lhsStride,"},
{"lineNum":"  109","line":"        const RhsScalar* _rhs, Index rhsIncr, ResScalar* _res, Index resIncr, const ResScalar& alpha)"},
{"lineNum":"  110","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  111","line":"    static const Index PanelWidth = EIGEN_TUNE_TRIANGULAR_PANEL_WIDTH;"},
{"lineNum":"  112","line":"    Index diagSize = (std::min)(_rows,_cols);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  113","line":"    Index rows = IsLower ? _rows : diagSize;"},
{"lineNum":"  114","line":"    Index cols = IsLower ? diagSize : _cols;"},
{"lineNum":"  115","line":""},
{"lineNum":"  116","line":"    typedef Map<const Matrix<LhsScalar,Dynamic,Dynamic,RowMajor>, 0, OuterStride<> > LhsMap;"},
{"lineNum":"  117","line":"    const LhsMap lhs(_lhs,rows,cols,OuterStride<>(lhsStride));"},
{"lineNum":"  118","line":"    typename conj_expr_if<ConjLhs,LhsMap>::type cjLhs(lhs);"},
{"lineNum":"  119","line":""},
{"lineNum":"  120","line":"    typedef Map<const Matrix<RhsScalar,Dynamic,1> > RhsMap;"},
{"lineNum":"  121","line":"    const RhsMap rhs(_rhs,cols);"},
{"lineNum":"  122","line":"    typename conj_expr_if<ConjRhs,RhsMap>::type cjRhs(rhs);"},
{"lineNum":"  123","line":""},
{"lineNum":"  124","line":"    typedef Map<Matrix<ResScalar,Dynamic,1>, 0, InnerStride<> > ResMap;"},
{"lineNum":"  125","line":"    ResMap res(_res,rows,InnerStride<>(resIncr));"},
{"lineNum":"  126","line":""},
{"lineNum":"  127","line":"    typedef const_blas_data_mapper<LhsScalar,Index,RowMajor> LhsMapper;"},
{"lineNum":"  128","line":"    typedef const_blas_data_mapper<RhsScalar,Index,RowMajor> RhsMapper;"},
{"lineNum":"  129","line":""},
{"lineNum":"  130","line":"    for (Index pi=0; pi<diagSize; pi+=PanelWidth)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  131","line":"    {"},
{"lineNum":"  132","line":"      Index actualPanelWidth = (std::min)(PanelWidth, diagSize-pi);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  133","line":"      for (Index k=0; k<actualPanelWidth; ++k)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  134","line":"      {"},
{"lineNum":"  135","line":"        Index i = pi + k;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  136","line":"        Index s = IsLower ? pi  : ((HasUnitDiag||HasZeroDiag) ? i+1 : i);"},
{"lineNum":"  137","line":"        Index r = IsLower ? k+1 : actualPanelWidth-k;"},
{"lineNum":"  138","line":"        if ((!(HasUnitDiag||HasZeroDiag)) || (--r)>0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  139","line":"          res.coeffRef(i) += alpha * (cjLhs.row(i).segment(s,r).cwiseProduct(cjRhs.segment(s,r).transpose())).sum();","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  140","line":"        if (HasUnitDiag)"},
{"lineNum":"  141","line":"          res.coeffRef(i) += alpha * cjRhs.coeff(i);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  142","line":"      }"},
{"lineNum":"  143","line":"      Index r = IsLower ? pi : cols - pi - actualPanelWidth;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  144","line":"      if (r>0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  145","line":"      {"},
{"lineNum":"  146","line":"        Index s = IsLower ? 0 : pi + actualPanelWidth;"},
{"lineNum":"  147","line":"        general_matrix_vector_product<Index,LhsScalar,LhsMapper,RowMajor,ConjLhs,RhsScalar,RhsMapper,ConjRhs,BuiltIn>::run(","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  148","line":"            actualPanelWidth, r,"},
{"lineNum":"  149","line":"            LhsMapper(&lhs.coeffRef(pi,s), lhsStride),"},
{"lineNum":"  150","line":"            RhsMapper(&rhs.coeffRef(s), rhsIncr),"},
{"lineNum":"  151","line":"            &res.coeffRef(pi), resIncr, alpha);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  152","line":"      }"},
{"lineNum":"  153","line":"    }"},
{"lineNum":"  154","line":"    if(IsLower && rows>diagSize)"},
{"lineNum":"  155","line":"    {"},
{"lineNum":"  156","line":"      general_matrix_vector_product<Index,LhsScalar,LhsMapper,RowMajor,ConjLhs,RhsScalar,RhsMapper,ConjRhs>::run("},
{"lineNum":"  157","line":"            rows-diagSize, cols,"},
{"lineNum":"  158","line":"            LhsMapper(&lhs.coeffRef(diagSize,0), lhsStride),"},
{"lineNum":"  159","line":"            RhsMapper(&rhs.coeffRef(0), rhsIncr),"},
{"lineNum":"  160","line":"            &res.coeffRef(diagSize), resIncr, alpha);"},
{"lineNum":"  161","line":"    }"},
{"lineNum":"  162","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  163","line":""},
{"lineNum":"  164","line":"/***************************************************************************"},
{"lineNum":"  165","line":"* Wrapper to product_triangular_vector"},
{"lineNum":"  166","line":"***************************************************************************/"},
{"lineNum":"  167","line":""},
{"lineNum":"  168","line":"template<int Mode,int StorageOrder>"},
{"lineNum":"  169","line":"struct trmv_selector;"},
{"lineNum":"  170","line":""},
{"lineNum":"  171","line":"} // end namespace internal"},
{"lineNum":"  172","line":""},
{"lineNum":"  173","line":"namespace internal {"},
{"lineNum":"  174","line":""},
{"lineNum":"  175","line":"template<int Mode, typename Lhs, typename Rhs>"},
{"lineNum":"  176","line":"struct triangular_product_impl<Mode,true,Lhs,false,Rhs,true>"},
{"lineNum":"  177","line":"{"},
{"lineNum":"  178","line":"  template<typename Dest> static void run(Dest& dst, const Lhs &lhs, const Rhs &rhs, const typename Dest::Scalar& alpha)"},
{"lineNum":"  179","line":"  {"},
{"lineNum":"  180","line":"    eigen_assert(dst.rows()==lhs.rows() && dst.cols()==rhs.cols());"},
{"lineNum":"  181","line":""},
{"lineNum":"  182","line":"    internal::trmv_selector<Mode,(int(internal::traits<Lhs>::Flags)&RowMajorBit) ? RowMajor : ColMajor>::run(lhs, rhs, dst, alpha);"},
{"lineNum":"  183","line":"  }"},
{"lineNum":"  184","line":"};"},
{"lineNum":"  185","line":""},
{"lineNum":"  186","line":"template<int Mode, typename Lhs, typename Rhs>"},
{"lineNum":"  187","line":"struct triangular_product_impl<Mode,false,Lhs,true,Rhs,false>"},
{"lineNum":"  188","line":"{"},
{"lineNum":"  189","line":"  template<typename Dest> static void run(Dest& dst, const Lhs &lhs, const Rhs &rhs, const typename Dest::Scalar& alpha)"},
{"lineNum":"  190","line":"  {"},
{"lineNum":"  191","line":"    eigen_assert(dst.rows()==lhs.rows() && dst.cols()==rhs.cols());"},
{"lineNum":"  192","line":""},
{"lineNum":"  193","line":"    Transpose<Dest> dstT(dst);"},
{"lineNum":"  194","line":"    internal::trmv_selector<(Mode & (UnitDiag|ZeroDiag)) | ((Mode & Lower) ? Upper : Lower),","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  195","line":"                            (int(internal::traits<Rhs>::Flags)&RowMajorBit) ? ColMajor : RowMajor>"},
{"lineNum":"  196","line":"            ::run(rhs.transpose(),lhs.transpose(), dstT, alpha);"},
{"lineNum":"  197","line":"  }"},
{"lineNum":"  198","line":"};"},
{"lineNum":"  199","line":""},
{"lineNum":"  200","line":"} // end namespace internal"},
{"lineNum":"  201","line":""},
{"lineNum":"  202","line":"namespace internal {"},
{"lineNum":"  203","line":""},
{"lineNum":"  204","line":"// TODO: find a way to factorize this piece of code with gemv_selector since the logic is exactly the same."},
{"lineNum":"  205","line":""},
{"lineNum":"  206","line":"template<int Mode> struct trmv_selector<Mode,ColMajor>"},
{"lineNum":"  207","line":"{"},
{"lineNum":"  208","line":"  template<typename Lhs, typename Rhs, typename Dest>"},
{"lineNum":"  209","line":"  static void run(const Lhs &lhs, const Rhs &rhs, Dest& dest, const typename Dest::Scalar& alpha)"},
{"lineNum":"  210","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  211","line":"    typedef typename Lhs::Scalar      LhsScalar;"},
{"lineNum":"  212","line":"    typedef typename Rhs::Scalar      RhsScalar;"},
{"lineNum":"  213","line":"    typedef typename Dest::Scalar     ResScalar;"},
{"lineNum":"  214","line":"    typedef typename Dest::RealScalar RealScalar;"},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"    typedef internal::blas_traits<Lhs> LhsBlasTraits;"},
{"lineNum":"  217","line":"    typedef typename LhsBlasTraits::DirectLinearAccessType ActualLhsType;"},
{"lineNum":"  218","line":"    typedef internal::blas_traits<Rhs> RhsBlasTraits;"},
{"lineNum":"  219","line":"    typedef typename RhsBlasTraits::DirectLinearAccessType ActualRhsType;"},
{"lineNum":"  220","line":""},
{"lineNum":"  221","line":"    typedef Map<Matrix<ResScalar,Dynamic,1>, EIGEN_PLAIN_ENUM_MIN(AlignedMax,internal::packet_traits<ResScalar>::size)> MappedDest;"},
{"lineNum":"  222","line":""},
{"lineNum":"  223","line":"    typename internal::add_const_on_value_type<ActualLhsType>::type actualLhs = LhsBlasTraits::extract(lhs);"},
{"lineNum":"  224","line":"    typename internal::add_const_on_value_type<ActualRhsType>::type actualRhs = RhsBlasTraits::extract(rhs);"},
{"lineNum":"  225","line":""},
{"lineNum":"  226","line":"    LhsScalar lhs_alpha = LhsBlasTraits::extractScalarFactor(lhs);"},
{"lineNum":"  227","line":"    RhsScalar rhs_alpha = RhsBlasTraits::extractScalarFactor(rhs);"},
{"lineNum":"  228","line":"    ResScalar actualAlpha = alpha * lhs_alpha * rhs_alpha;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  229","line":""},
{"lineNum":"  230","line":"    enum {"},
{"lineNum":"  231","line":"      // FIXME find a way to allow an inner stride on the result if packet_traits<Scalar>::size==1"},
{"lineNum":"  232","line":"      // on, the other hand it is good for the cache to pack the vector anyways..."},
{"lineNum":"  233","line":"      EvalToDestAtCompileTime = Dest::InnerStrideAtCompileTime==1,"},
{"lineNum":"  234","line":"      ComplexByReal = (NumTraits<LhsScalar>::IsComplex) && (!NumTraits<RhsScalar>::IsComplex),"},
{"lineNum":"  235","line":"      MightCannotUseDest = (Dest::InnerStrideAtCompileTime!=1) || ComplexByReal"},
{"lineNum":"  236","line":"    };"},
{"lineNum":"  237","line":""},
{"lineNum":"  238","line":"    gemv_static_vector_if<ResScalar,Dest::SizeAtCompileTime,Dest::MaxSizeAtCompileTime,MightCannotUseDest> static_dest;"},
{"lineNum":"  239","line":""},
{"lineNum":"  240","line":"    bool alphaIsCompatible = (!ComplexByReal) || (numext::imag(actualAlpha)==RealScalar(0));"},
{"lineNum":"  241","line":"    bool evalToDest = EvalToDestAtCompileTime && alphaIsCompatible;"},
{"lineNum":"  242","line":""},
{"lineNum":"  243","line":"    RhsScalar compatibleAlpha = get_factor<ResScalar,RhsScalar>::run(actualAlpha);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"    ei_declare_aligned_stack_constructed_variable(ResScalar,actualDestPtr,dest.size(),","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  246","line":"                                                  evalToDest ? dest.data() : static_dest.data());"},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"    if(!evalToDest)"},
{"lineNum":"  249","line":"    {"},
{"lineNum":"  250","line":"      #ifdef EIGEN_DENSE_STORAGE_CTOR_PLUGIN"},
{"lineNum":"  251","line":"      Index size = dest.size();"},
{"lineNum":"  252","line":"      EIGEN_DENSE_STORAGE_CTOR_PLUGIN"},
{"lineNum":"  253","line":"      #endif"},
{"lineNum":"  254","line":"      if(!alphaIsCompatible)"},
{"lineNum":"  255","line":"      {"},
{"lineNum":"  256","line":"        MappedDest(actualDestPtr, dest.size()).setZero();"},
{"lineNum":"  257","line":"        compatibleAlpha = RhsScalar(1);"},
{"lineNum":"  258","line":"      }"},
{"lineNum":"  259","line":"      else"},
{"lineNum":"  260","line":"        MappedDest(actualDestPtr, dest.size()) = dest;"},
{"lineNum":"  261","line":"    }"},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"    internal::triangular_matrix_vector_product","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  264","line":"      <Index,Mode,"},
{"lineNum":"  265","line":"       LhsScalar, LhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  266","line":"       RhsScalar, RhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  267","line":"       ColMajor>"},
{"lineNum":"  268","line":"      ::run(actualLhs.rows(),actualLhs.cols(),"},
{"lineNum":"  269","line":"            actualLhs.data(),actualLhs.outerStride(),"},
{"lineNum":"  270","line":"            actualRhs.data(),actualRhs.innerStride(),"},
{"lineNum":"  271","line":"            actualDestPtr,1,compatibleAlpha);"},
{"lineNum":"  272","line":""},
{"lineNum":"  273","line":"    if (!evalToDest)"},
{"lineNum":"  274","line":"    {"},
{"lineNum":"  275","line":"      if(!alphaIsCompatible)"},
{"lineNum":"  276","line":"        dest += actualAlpha * MappedDest(actualDestPtr, dest.size());"},
{"lineNum":"  277","line":"      else"},
{"lineNum":"  278","line":"        dest = MappedDest(actualDestPtr, dest.size());"},
{"lineNum":"  279","line":"    }"},
{"lineNum":"  280","line":""},
{"lineNum":"  281","line":"    if ( ((Mode&UnitDiag)==UnitDiag) && (lhs_alpha!=LhsScalar(1)) )"},
{"lineNum":"  282","line":"    {"},
{"lineNum":"  283","line":"      Index diagSize = (std::min)(lhs.rows(),lhs.cols());"},
{"lineNum":"  284","line":"      dest.head(diagSize) -= (lhs_alpha-LhsScalar(1))*rhs.head(diagSize);"},
{"lineNum":"  285","line":"    }"},
{"lineNum":"  286","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  287","line":"};"},
{"lineNum":"  288","line":""},
{"lineNum":"  289","line":"template<int Mode> struct trmv_selector<Mode,RowMajor>"},
{"lineNum":"  290","line":"{"},
{"lineNum":"  291","line":"  template<typename Lhs, typename Rhs, typename Dest>"},
{"lineNum":"  292","line":"  static void run(const Lhs &lhs, const Rhs &rhs, Dest& dest, const typename Dest::Scalar& alpha)"},
{"lineNum":"  293","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  294","line":"    typedef typename Lhs::Scalar      LhsScalar;"},
{"lineNum":"  295","line":"    typedef typename Rhs::Scalar      RhsScalar;"},
{"lineNum":"  296","line":"    typedef typename Dest::Scalar     ResScalar;"},
{"lineNum":"  297","line":""},
{"lineNum":"  298","line":"    typedef internal::blas_traits<Lhs> LhsBlasTraits;"},
{"lineNum":"  299","line":"    typedef typename LhsBlasTraits::DirectLinearAccessType ActualLhsType;"},
{"lineNum":"  300","line":"    typedef internal::blas_traits<Rhs> RhsBlasTraits;"},
{"lineNum":"  301","line":"    typedef typename RhsBlasTraits::DirectLinearAccessType ActualRhsType;"},
{"lineNum":"  302","line":"    typedef typename internal::remove_all<ActualRhsType>::type ActualRhsTypeCleaned;"},
{"lineNum":"  303","line":""},
{"lineNum":"  304","line":"    typename add_const<ActualLhsType>::type actualLhs = LhsBlasTraits::extract(lhs);"},
{"lineNum":"  305","line":"    typename add_const<ActualRhsType>::type actualRhs = RhsBlasTraits::extract(rhs);"},
{"lineNum":"  306","line":""},
{"lineNum":"  307","line":"    LhsScalar lhs_alpha = LhsBlasTraits::extractScalarFactor(lhs);"},
{"lineNum":"  308","line":"    RhsScalar rhs_alpha = RhsBlasTraits::extractScalarFactor(rhs);"},
{"lineNum":"  309","line":"    ResScalar actualAlpha = alpha * lhs_alpha * rhs_alpha;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  310","line":""},
{"lineNum":"  311","line":"    enum {"},
{"lineNum":"  312","line":"      DirectlyUseRhs = ActualRhsTypeCleaned::InnerStrideAtCompileTime==1"},
{"lineNum":"  313","line":"    };"},
{"lineNum":"  314","line":""},
{"lineNum":"  315","line":"    gemv_static_vector_if<RhsScalar,ActualRhsTypeCleaned::SizeAtCompileTime,ActualRhsTypeCleaned::MaxSizeAtCompileTime,!DirectlyUseRhs> static_rhs;"},
{"lineNum":"  316","line":""},
{"lineNum":"  317","line":"    ei_declare_aligned_stack_constructed_variable(RhsScalar,actualRhsPtr,actualRhs.size(),","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  318","line":"        DirectlyUseRhs ? const_cast<RhsScalar*>(actualRhs.data()) : static_rhs.data());"},
{"lineNum":"  319","line":""},
{"lineNum":"  320","line":"    if(!DirectlyUseRhs)"},
{"lineNum":"  321","line":"    {"},
{"lineNum":"  322","line":"      #ifdef EIGEN_DENSE_STORAGE_CTOR_PLUGIN"},
{"lineNum":"  323","line":"      Index size = actualRhs.size();"},
{"lineNum":"  324","line":"      EIGEN_DENSE_STORAGE_CTOR_PLUGIN"},
{"lineNum":"  325","line":"      #endif"},
{"lineNum":"  326","line":"      Map<typename ActualRhsTypeCleaned::PlainObject>(actualRhsPtr, actualRhs.size()) = actualRhs;"},
{"lineNum":"  327","line":"    }"},
{"lineNum":"  328","line":""},
{"lineNum":"  329","line":"    internal::triangular_matrix_vector_product","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  330","line":"      <Index,Mode,"},
{"lineNum":"  331","line":"       LhsScalar, LhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  332","line":"       RhsScalar, RhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  333","line":"       RowMajor>"},
{"lineNum":"  334","line":"      ::run(actualLhs.rows(),actualLhs.cols(),"},
{"lineNum":"  335","line":"            actualLhs.data(),actualLhs.outerStride(),"},
{"lineNum":"  336","line":"            actualRhsPtr,1,"},
{"lineNum":"  337","line":"            dest.data(),dest.innerStride(),"},
{"lineNum":"  338","line":"            actualAlpha);"},
{"lineNum":"  339","line":""},
{"lineNum":"  340","line":"    if ( ((Mode&UnitDiag)==UnitDiag) && (lhs_alpha!=LhsScalar(1)) )"},
{"lineNum":"  341","line":"    {"},
{"lineNum":"  342","line":"      Index diagSize = (std::min)(lhs.rows(),lhs.cols());"},
{"lineNum":"  343","line":"      dest.head(diagSize) -= (lhs_alpha-LhsScalar(1))*rhs.head(diagSize);"},
{"lineNum":"  344","line":"    }"},
{"lineNum":"  345","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  346","line":"};"},
{"lineNum":"  347","line":""},
{"lineNum":"  348","line":"} // end namespace internal"},
{"lineNum":"  349","line":""},
{"lineNum":"  350","line":"} // end namespace Eigen"},
{"lineNum":"  351","line":""},
{"lineNum":"  352","line":"#endif // EIGEN_TRIANGULARMATRIXVECTOR_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 44, "covered" : 0,};
var merged_data = [];
