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
{"lineNum":"   10","line":"#ifndef EIGEN_TRIANGULAR_MATRIX_MATRIX_H"},
{"lineNum":"   11","line":"#define EIGEN_TRIANGULAR_MATRIX_MATRIX_H"},
{"lineNum":"   12","line":""},
{"lineNum":"   13","line":"#include \"../InternalHeaderCheck.h\""},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"namespace Eigen {"},
{"lineNum":"   16","line":""},
{"lineNum":"   17","line":"namespace internal {"},
{"lineNum":"   18","line":""},
{"lineNum":"   19","line":"// template<typename Scalar, int mr, int StorageOrder, bool Conjugate, int Mode>"},
{"lineNum":"   20","line":"// struct gemm_pack_lhs_triangular"},
{"lineNum":"   21","line":"// {"},
{"lineNum":"   22","line":"//   Matrix<Scalar,mr,mr,"},
{"lineNum":"   23","line":"//   void operator()(Scalar* blockA, const EIGEN_RESTRICT Scalar* _lhs, int lhsStride, int depth, int rows)"},
{"lineNum":"   24","line":"//   {"},
{"lineNum":"   25","line":"//     conj_if<NumTraits<Scalar>::IsComplex && Conjugate> cj;"},
{"lineNum":"   26","line":"//     const_blas_data_mapper<Scalar, StorageOrder> lhs(_lhs,lhsStride);"},
{"lineNum":"   27","line":"//     int count = 0;"},
{"lineNum":"   28","line":"//     const int peeled_mc = (rows/mr)*mr;"},
{"lineNum":"   29","line":"//     for(int i=0; i<peeled_mc; i+=mr)"},
{"lineNum":"   30","line":"//     {"},
{"lineNum":"   31","line":"//       for(int k=0; k<depth; k++)"},
{"lineNum":"   32","line":"//         for(int w=0; w<mr; w++)"},
{"lineNum":"   33","line":"//           blockA[count++] = cj(lhs(i+w, k));"},
{"lineNum":"   34","line":"//     }"},
{"lineNum":"   35","line":"//     for(int i=peeled_mc; i<rows; i++)"},
{"lineNum":"   36","line":"//     {"},
{"lineNum":"   37","line":"//       for(int k=0; k<depth; k++)"},
{"lineNum":"   38","line":"//         blockA[count++] = cj(lhs(i, k));"},
{"lineNum":"   39","line":"//     }"},
{"lineNum":"   40","line":"//   }"},
{"lineNum":"   41","line":"// };"},
{"lineNum":"   42","line":""},
{"lineNum":"   43","line":"/* Optimized triangular matrix * matrix (_TRMM++) product built on top of"},
{"lineNum":"   44","line":" * the general matrix matrix product."},
{"lineNum":"   45","line":" */"},
{"lineNum":"   46","line":"template <typename Scalar, typename Index,"},
{"lineNum":"   47","line":"          int Mode, bool LhsIsTriangular,"},
{"lineNum":"   48","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"   49","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"   50","line":"          int ResStorageOrder, int ResInnerStride,"},
{"lineNum":"   51","line":"          int Version = Specialized>"},
{"lineNum":"   52","line":"struct product_triangular_matrix_matrix;"},
{"lineNum":"   53","line":""},
{"lineNum":"   54","line":"template <typename Scalar, typename Index,"},
{"lineNum":"   55","line":"          int Mode, bool LhsIsTriangular,"},
{"lineNum":"   56","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"   57","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"   58","line":"          int ResInnerStride, int Version>"},
{"lineNum":"   59","line":"struct product_triangular_matrix_matrix<Scalar,Index,Mode,LhsIsTriangular,"},
{"lineNum":"   60","line":"                                           LhsStorageOrder,ConjugateLhs,"},
{"lineNum":"   61","line":"                                           RhsStorageOrder,ConjugateRhs,RowMajor,ResInnerStride,Version>"},
{"lineNum":"   62","line":"{"},
{"lineNum":"   63","line":"  static EIGEN_STRONG_INLINE void run("},
{"lineNum":"   64","line":"    Index rows, Index cols, Index depth,"},
{"lineNum":"   65","line":"    const Scalar* lhs, Index lhsStride,"},
{"lineNum":"   66","line":"    const Scalar* rhs, Index rhsStride,"},
{"lineNum":"   67","line":"    Scalar* res,       Index resIncr, Index resStride,"},
{"lineNum":"   68","line":"    const Scalar& alpha, level3_blocking<Scalar,Scalar>& blocking)"},
{"lineNum":"   69","line":"  {"},
{"lineNum":"   70","line":"    product_triangular_matrix_matrix<Scalar, Index,","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"   71","line":"      (Mode&(UnitDiag|ZeroDiag)) | ((Mode&Upper) ? Lower : Upper),"},
{"lineNum":"   72","line":"      (!LhsIsTriangular),"},
{"lineNum":"   73","line":"      RhsStorageOrder==RowMajor ? ColMajor : RowMajor,"},
{"lineNum":"   74","line":"      ConjugateRhs,"},
{"lineNum":"   75","line":"      LhsStorageOrder==RowMajor ? ColMajor : RowMajor,"},
{"lineNum":"   76","line":"      ConjugateLhs,"},
{"lineNum":"   77","line":"      ColMajor, ResInnerStride>"},
{"lineNum":"   78","line":"      ::run(cols, rows, depth, rhs, rhsStride, lhs, lhsStride, res, resIncr, resStride, alpha, blocking);"},
{"lineNum":"   79","line":"  }"},
{"lineNum":"   80","line":"};"},
{"lineNum":"   81","line":""},
{"lineNum":"   82","line":"// implements col-major += alpha * op(triangular) * op(general)"},
{"lineNum":"   83","line":"template <typename Scalar, typename Index, int Mode,"},
{"lineNum":"   84","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"   85","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"   86","line":"          int ResInnerStride, int Version>"},
{"lineNum":"   87","line":"struct product_triangular_matrix_matrix<Scalar,Index,Mode,true,"},
{"lineNum":"   88","line":"                                           LhsStorageOrder,ConjugateLhs,"},
{"lineNum":"   89","line":"                                           RhsStorageOrder,ConjugateRhs,ColMajor,ResInnerStride,Version>"},
{"lineNum":"   90","line":"{"},
{"lineNum":"   91","line":""},
{"lineNum":"   92","line":"  typedef gebp_traits<Scalar,Scalar> Traits;"},
{"lineNum":"   93","line":"  enum {"},
{"lineNum":"   94","line":"    SmallPanelWidth   = 2 * EIGEN_PLAIN_ENUM_MAX(Traits::mr,Traits::nr),"},
{"lineNum":"   95","line":"    IsLower = (Mode&Lower) == Lower,"},
{"lineNum":"   96","line":"    SetDiag = (Mode&(ZeroDiag|UnitDiag)) ? 0 : 1"},
{"lineNum":"   97","line":"  };"},
{"lineNum":"   98","line":""},
{"lineNum":"   99","line":"  static EIGEN_DONT_INLINE void run("},
{"lineNum":"  100","line":"    Index _rows, Index _cols, Index _depth,"},
{"lineNum":"  101","line":"    const Scalar* _lhs, Index lhsStride,"},
{"lineNum":"  102","line":"    const Scalar* _rhs, Index rhsStride,"},
{"lineNum":"  103","line":"    Scalar* res,        Index resIncr, Index resStride,"},
{"lineNum":"  104","line":"    const Scalar& alpha, level3_blocking<Scalar,Scalar>& blocking);"},
{"lineNum":"  105","line":"};"},
{"lineNum":"  106","line":""},
{"lineNum":"  107","line":"template <typename Scalar, typename Index, int Mode,"},
{"lineNum":"  108","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"  109","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"  110","line":"          int ResInnerStride, int Version>"},
{"lineNum":"  111","line":"EIGEN_DONT_INLINE void product_triangular_matrix_matrix<Scalar,Index,Mode,true,"},
{"lineNum":"  112","line":"                                                        LhsStorageOrder,ConjugateLhs,"},
{"lineNum":"  113","line":"                                                        RhsStorageOrder,ConjugateRhs,ColMajor,ResInnerStride,Version>::run("},
{"lineNum":"  114","line":"    Index _rows, Index _cols, Index _depth,"},
{"lineNum":"  115","line":"    const Scalar* _lhs, Index lhsStride,"},
{"lineNum":"  116","line":"    const Scalar* _rhs, Index rhsStride,"},
{"lineNum":"  117","line":"    Scalar* _res,       Index resIncr, Index resStride,"},
{"lineNum":"  118","line":"    const Scalar& alpha, level3_blocking<Scalar,Scalar>& blocking)"},
{"lineNum":"  119","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  120","line":"    // strip zeros"},
{"lineNum":"  121","line":"    Index diagSize  = (std::min)(_rows,_depth);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  122","line":"    Index rows      = IsLower ? _rows : diagSize;"},
{"lineNum":"  123","line":"    Index depth     = IsLower ? diagSize : _depth;"},
{"lineNum":"  124","line":"    Index cols      = _cols;"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"    typedef const_blas_data_mapper<Scalar, Index, LhsStorageOrder> LhsMapper;"},
{"lineNum":"  127","line":"    typedef const_blas_data_mapper<Scalar, Index, RhsStorageOrder> RhsMapper;"},
{"lineNum":"  128","line":"    typedef blas_data_mapper<typename Traits::ResScalar, Index, ColMajor, Unaligned, ResInnerStride> ResMapper;"},
{"lineNum":"  129","line":"    LhsMapper lhs(_lhs,lhsStride);"},
{"lineNum":"  130","line":"    RhsMapper rhs(_rhs,rhsStride);"},
{"lineNum":"  131","line":"    ResMapper res(_res, resStride, resIncr);"},
{"lineNum":"  132","line":""},
{"lineNum":"  133","line":"    Index kc = blocking.kc();                   // cache block size along the K direction"},
{"lineNum":"  134","line":"    Index mc = (std::min)(rows,blocking.mc());  // cache block size along the M direction","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  135","line":"    // The small panel size must not be larger than blocking size."},
{"lineNum":"  136","line":"    // Usually this should never be the case because SmallPanelWidth^2 is very small"},
{"lineNum":"  137","line":"    // compared to L2 cache size, but let\'s be safe:"},
{"lineNum":"  138","line":"    Index panelWidth = (std::min)(Index(SmallPanelWidth),(std::min)(kc,mc));","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  139","line":""},
{"lineNum":"  140","line":"    std::size_t sizeA = kc*mc;","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  141","line":"    std::size_t sizeB = kc*cols;"},
{"lineNum":"  142","line":""},
{"lineNum":"  143","line":"    ei_declare_aligned_stack_constructed_variable(Scalar, blockA, sizeA, blocking.blockA());","class":"lineNoCov","hits":"0","possible_hits":"21",},
{"lineNum":"  144","line":"    ei_declare_aligned_stack_constructed_variable(Scalar, blockB, sizeB, blocking.blockB());","class":"lineNoCov","hits":"0","possible_hits":"21",},
{"lineNum":"  145","line":""},
{"lineNum":"  146","line":"    // To work around an \"error: member reference base type \'Matrix<...>"},
{"lineNum":"  147","line":"    // (Eigen::internal::constructor_without_unaligned_array_assert (*)())\' is"},
{"lineNum":"  148","line":"    // not a structure or union\" compilation error in nvcc (tested V8.0.61),"},
{"lineNum":"  149","line":"    // create a dummy internal::constructor_without_unaligned_array_assert"},
{"lineNum":"  150","line":"    // object to pass to the Matrix constructor."},
{"lineNum":"  151","line":"    internal::constructor_without_unaligned_array_assert a;"},
{"lineNum":"  152","line":"    Matrix<Scalar,SmallPanelWidth,SmallPanelWidth,LhsStorageOrder> triangularBuffer(a);"},
{"lineNum":"  153","line":"    triangularBuffer.setZero();"},
{"lineNum":"  154","line":"    if((Mode&ZeroDiag)==ZeroDiag)"},
{"lineNum":"  155","line":"      triangularBuffer.diagonal().setZero();"},
{"lineNum":"  156","line":"    else"},
{"lineNum":"  157","line":"      triangularBuffer.diagonal().setOnes();"},
{"lineNum":"  158","line":""},
{"lineNum":"  159","line":"    gebp_kernel<Scalar, Scalar, Index, ResMapper, Traits::mr, Traits::nr, ConjugateLhs, ConjugateRhs> gebp_kernel;"},
{"lineNum":"  160","line":"    gemm_pack_lhs<Scalar, Index, LhsMapper, Traits::mr, Traits::LhsProgress, typename Traits::LhsPacket4Packing, LhsStorageOrder> pack_lhs;"},
{"lineNum":"  161","line":"    gemm_pack_rhs<Scalar, Index, RhsMapper, Traits::nr,RhsStorageOrder> pack_rhs;"},
{"lineNum":"  162","line":""},
{"lineNum":"  163","line":"    for(Index k2=IsLower ? depth : 0;","class":"lineNoCov","hits":"0","possible_hits":"23",},
{"lineNum":"  164","line":"        IsLower ? k2>0 : k2<depth;","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  165","line":"        IsLower ? k2-=kc : k2+=kc)","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  166","line":"    {"},
{"lineNum":"  167","line":"      Index actual_kc = (std::min)(IsLower ? k2 : depth-k2, kc);","class":"lineNoCov","hits":"0","possible_hits":"11",},
{"lineNum":"  168","line":"      Index actual_k2 = IsLower ? k2-actual_kc : k2;","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  169","line":""},
{"lineNum":"  170","line":"      // align blocks with the end of the triangular part for trapezoidal lhs"},
{"lineNum":"  171","line":"      if((!IsLower)&&(k2<rows)&&(k2+actual_kc>rows))","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"  172","line":"      {"},
{"lineNum":"  173","line":"        actual_kc = rows-k2;"},
{"lineNum":"  174","line":"        k2 = k2+actual_kc-kc;"},
{"lineNum":"  175","line":"      }"},
{"lineNum":"  176","line":""},
{"lineNum":"  177","line":"      pack_rhs(blockB, rhs.getSubMapper(actual_k2,0), actual_kc, cols);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  178","line":""},
{"lineNum":"  179","line":"      // the selected lhs\'s panel has to be split in three different parts:"},
{"lineNum":"  180","line":"      //  1 - the part which is zero => skip it"},
{"lineNum":"  181","line":"      //  2 - the diagonal block => special kernel"},
{"lineNum":"  182","line":"      //  3 - the dense panel below (lower case) or above (upper case) the diagonal block => GEPP"},
{"lineNum":"  183","line":""},
{"lineNum":"  184","line":"      // the block diagonal, if any:"},
{"lineNum":"  185","line":"      if(IsLower || actual_k2<rows)","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  186","line":"      {"},
{"lineNum":"  187","line":"        // for each small vertical panels of lhs"},
{"lineNum":"  188","line":"        for (Index k1=0; k1<actual_kc; k1+=panelWidth)","class":"lineNoCov","hits":"0","possible_hits":"16",},
{"lineNum":"  189","line":"        {"},
{"lineNum":"  190","line":"          Index actualPanelWidth = std::min<Index>(actual_kc-k1, panelWidth);","class":"lineNoCov","hits":"0","possible_hits":"16",},
{"lineNum":"  191","line":"          Index lengthTarget = IsLower ? actual_kc-k1-actualPanelWidth : k1;"},
{"lineNum":"  192","line":"          Index startBlock   = actual_k2+k1;","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  193","line":"          Index blockBOffset = k1;"},
{"lineNum":"  194","line":""},
{"lineNum":"  195","line":"          // => GEBP with the micro triangular block"},
{"lineNum":"  196","line":"          // The trick is to pack this micro block while filling the opposite triangular part with zeros."},
{"lineNum":"  197","line":"          // To this end we do an extra triangular copy to a small temporary buffer"},
{"lineNum":"  198","line":"          for (Index k=0;k<actualPanelWidth;++k)","class":"lineNoCov","hits":"0","possible_hits":"9",},
{"lineNum":"  199","line":"          {"},
{"lineNum":"  200","line":"            if (SetDiag)"},
{"lineNum":"  201","line":"              triangularBuffer.coeffRef(k,k) = lhs(startBlock+k,startBlock+k);","class":"lineNoCov","hits":"0","possible_hits":"10",},
{"lineNum":"  202","line":"            for (Index i=IsLower ? k+1 : 0; IsLower ? i<actualPanelWidth : i<k; ++i)","class":"lineNoCov","hits":"0","possible_hits":"12",},
{"lineNum":"  203","line":"              triangularBuffer.coeffRef(i,k) = lhs(startBlock+i,startBlock+k);","class":"lineNoCov","hits":"0","possible_hits":"19",},
{"lineNum":"  204","line":"          }"},
{"lineNum":"  205","line":"          pack_lhs(blockA, LhsMapper(triangularBuffer.data(), triangularBuffer.outerStride()), actualPanelWidth, actualPanelWidth);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  206","line":""},
{"lineNum":"  207","line":"          gebp_kernel(res.getSubMapper(startBlock, 0), blockA, blockB,","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  208","line":"                      actualPanelWidth, actualPanelWidth, cols, alpha,","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  209","line":"                      actualPanelWidth, actual_kc, 0, blockBOffset);"},
{"lineNum":"  210","line":""},
{"lineNum":"  211","line":"          // GEBP with remaining micro panel"},
{"lineNum":"  212","line":"          if (lengthTarget>0)","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  213","line":"          {"},
{"lineNum":"  214","line":"            Index startTarget  = IsLower ? actual_k2+k1+actualPanelWidth : actual_k2;","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"            pack_lhs(blockA, lhs.getSubMapper(startTarget,startBlock), actualPanelWidth, lengthTarget);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  217","line":""},
{"lineNum":"  218","line":"            gebp_kernel(res.getSubMapper(startTarget, 0), blockA, blockB,","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  219","line":"                        lengthTarget, actualPanelWidth, cols, alpha,","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  220","line":"                        actualPanelWidth, actual_kc, 0, blockBOffset);"},
{"lineNum":"  221","line":"          }"},
{"lineNum":"  222","line":"        }"},
{"lineNum":"  223","line":"      }"},
{"lineNum":"  224","line":"      // the part below (lower case) or above (upper case) the diagonal => GEPP"},
{"lineNum":"  225","line":"      {"},
{"lineNum":"  226","line":"        Index start = IsLower ? k2 : 0;"},
{"lineNum":"  227","line":"        Index end   = IsLower ? rows : (std::min)(actual_k2,rows);","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  228","line":"        for(Index i2=start; i2<end; i2+=mc)","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  229","line":"        {"},
{"lineNum":"  230","line":"          const Index actual_mc = (std::min)(i2+mc,end)-i2;","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  231","line":"          gemm_pack_lhs<Scalar, Index, LhsMapper, Traits::mr,Traits::LhsProgress, typename Traits::LhsPacket4Packing, LhsStorageOrder,false>()","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  232","line":"            (blockA, lhs.getSubMapper(i2, actual_k2), actual_kc, actual_mc);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  233","line":""},
{"lineNum":"  234","line":"          gebp_kernel(res.getSubMapper(i2, 0), blockA, blockB, actual_mc,","class":"lineNoCov","hits":"0","possible_hits":"14",},
{"lineNum":"  235","line":"                      actual_kc, cols, alpha, -1, -1, 0, 0);","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  236","line":"        }"},
{"lineNum":"  237","line":"      }"},
{"lineNum":"  238","line":"    }"},
{"lineNum":"  239","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  240","line":""},
{"lineNum":"  241","line":"// implements col-major += alpha * op(general) * op(triangular)"},
{"lineNum":"  242","line":"template <typename Scalar, typename Index, int Mode,"},
{"lineNum":"  243","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"  244","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"  245","line":"          int ResInnerStride, int Version>"},
{"lineNum":"  246","line":"struct product_triangular_matrix_matrix<Scalar,Index,Mode,false,"},
{"lineNum":"  247","line":"                                        LhsStorageOrder,ConjugateLhs,"},
{"lineNum":"  248","line":"                                        RhsStorageOrder,ConjugateRhs,ColMajor,ResInnerStride,Version>"},
{"lineNum":"  249","line":"{"},
{"lineNum":"  250","line":"  typedef gebp_traits<Scalar,Scalar> Traits;"},
{"lineNum":"  251","line":"  enum {"},
{"lineNum":"  252","line":"    SmallPanelWidth   = EIGEN_PLAIN_ENUM_MAX(Traits::mr,Traits::nr),"},
{"lineNum":"  253","line":"    IsLower = (Mode&Lower) == Lower,"},
{"lineNum":"  254","line":"    SetDiag = (Mode&(ZeroDiag|UnitDiag)) ? 0 : 1"},
{"lineNum":"  255","line":"  };"},
{"lineNum":"  256","line":""},
{"lineNum":"  257","line":"  static EIGEN_DONT_INLINE void run("},
{"lineNum":"  258","line":"    Index _rows, Index _cols, Index _depth,"},
{"lineNum":"  259","line":"    const Scalar* _lhs, Index lhsStride,"},
{"lineNum":"  260","line":"    const Scalar* _rhs, Index rhsStride,"},
{"lineNum":"  261","line":"    Scalar* res,        Index resIncr, Index resStride,"},
{"lineNum":"  262","line":"    const Scalar& alpha, level3_blocking<Scalar,Scalar>& blocking);"},
{"lineNum":"  263","line":"};"},
{"lineNum":"  264","line":""},
{"lineNum":"  265","line":"template <typename Scalar, typename Index, int Mode,"},
{"lineNum":"  266","line":"          int LhsStorageOrder, bool ConjugateLhs,"},
{"lineNum":"  267","line":"          int RhsStorageOrder, bool ConjugateRhs,"},
{"lineNum":"  268","line":"          int ResInnerStride, int Version>"},
{"lineNum":"  269","line":"EIGEN_DONT_INLINE void product_triangular_matrix_matrix<Scalar,Index,Mode,false,"},
{"lineNum":"  270","line":"                                                        LhsStorageOrder,ConjugateLhs,"},
{"lineNum":"  271","line":"                                                        RhsStorageOrder,ConjugateRhs,ColMajor,ResInnerStride,Version>::run("},
{"lineNum":"  272","line":"    Index _rows, Index _cols, Index _depth,"},
{"lineNum":"  273","line":"    const Scalar* _lhs, Index lhsStride,"},
{"lineNum":"  274","line":"    const Scalar* _rhs, Index rhsStride,"},
{"lineNum":"  275","line":"    Scalar* _res,       Index resIncr, Index resStride,"},
{"lineNum":"  276","line":"    const Scalar& alpha, level3_blocking<Scalar,Scalar>& blocking)"},
{"lineNum":"  277","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  278","line":"    const Index PacketBytes = packet_traits<Scalar>::size*sizeof(Scalar);"},
{"lineNum":"  279","line":"    // strip zeros"},
{"lineNum":"  280","line":"    Index diagSize  = (std::min)(_cols,_depth);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  281","line":"    Index rows      = _rows;"},
{"lineNum":"  282","line":"    Index depth     = IsLower ? _depth : diagSize;"},
{"lineNum":"  283","line":"    Index cols      = IsLower ? diagSize : _cols;"},
{"lineNum":"  284","line":""},
{"lineNum":"  285","line":"    typedef const_blas_data_mapper<Scalar, Index, LhsStorageOrder> LhsMapper;"},
{"lineNum":"  286","line":"    typedef const_blas_data_mapper<Scalar, Index, RhsStorageOrder> RhsMapper;"},
{"lineNum":"  287","line":"    typedef blas_data_mapper<typename Traits::ResScalar, Index, ColMajor, Unaligned, ResInnerStride> ResMapper;"},
{"lineNum":"  288","line":"    LhsMapper lhs(_lhs,lhsStride);"},
{"lineNum":"  289","line":"    RhsMapper rhs(_rhs,rhsStride);"},
{"lineNum":"  290","line":"    ResMapper res(_res, resStride, resIncr);"},
{"lineNum":"  291","line":""},
{"lineNum":"  292","line":"    Index kc = blocking.kc();                   // cache block size along the K direction"},
{"lineNum":"  293","line":"    Index mc = (std::min)(rows,blocking.mc());  // cache block size along the M direction","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  294","line":""},
{"lineNum":"  295","line":"    std::size_t sizeA = kc*mc;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  296","line":"    std::size_t sizeB = kc*cols+EIGEN_MAX_ALIGN_BYTES/sizeof(Scalar);"},
{"lineNum":"  297","line":""},
{"lineNum":"  298","line":"    ei_declare_aligned_stack_constructed_variable(Scalar, blockA, sizeA, blocking.blockA());","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  299","line":"    ei_declare_aligned_stack_constructed_variable(Scalar, blockB, sizeB, blocking.blockB());","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  300","line":""},
{"lineNum":"  301","line":"    internal::constructor_without_unaligned_array_assert a;"},
{"lineNum":"  302","line":"    Matrix<Scalar,SmallPanelWidth,SmallPanelWidth,RhsStorageOrder> triangularBuffer(a);"},
{"lineNum":"  303","line":"    triangularBuffer.setZero();"},
{"lineNum":"  304","line":"    if((Mode&ZeroDiag)==ZeroDiag)"},
{"lineNum":"  305","line":"      triangularBuffer.diagonal().setZero();"},
{"lineNum":"  306","line":"    else"},
{"lineNum":"  307","line":"      triangularBuffer.diagonal().setOnes();"},
{"lineNum":"  308","line":""},
{"lineNum":"  309","line":"    gebp_kernel<Scalar, Scalar, Index, ResMapper, Traits::mr, Traits::nr, ConjugateLhs, ConjugateRhs> gebp_kernel;"},
{"lineNum":"  310","line":"    gemm_pack_lhs<Scalar, Index, LhsMapper, Traits::mr, Traits::LhsProgress, typename Traits::LhsPacket4Packing, LhsStorageOrder> pack_lhs;"},
{"lineNum":"  311","line":"    gemm_pack_rhs<Scalar, Index, RhsMapper, Traits::nr,RhsStorageOrder> pack_rhs;"},
{"lineNum":"  312","line":"    gemm_pack_rhs<Scalar, Index, RhsMapper, Traits::nr,RhsStorageOrder,false,true> pack_rhs_panel;"},
{"lineNum":"  313","line":""},
{"lineNum":"  314","line":"    for(Index k2=IsLower ? 0 : depth;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  315","line":"        IsLower ? k2<depth  : k2>0;","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  316","line":"        IsLower ? k2+=kc   : k2-=kc)"},
{"lineNum":"  317","line":"    {"},
{"lineNum":"  318","line":"      Index actual_kc = (std::min)(IsLower ? depth-k2 : k2, kc);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  319","line":"      Index actual_k2 = IsLower ? k2 : k2-actual_kc;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  320","line":""},
{"lineNum":"  321","line":"      // align blocks with the end of the triangular part for trapezoidal rhs"},
{"lineNum":"  322","line":"      if(IsLower && (k2<cols) && (actual_k2+actual_kc>cols))"},
{"lineNum":"  323","line":"      {"},
{"lineNum":"  324","line":"        actual_kc = cols-k2;"},
{"lineNum":"  325","line":"        k2 = actual_k2 + actual_kc - kc;"},
{"lineNum":"  326","line":"      }"},
{"lineNum":"  327","line":""},
{"lineNum":"  328","line":"      // remaining size"},
{"lineNum":"  329","line":"      Index rs = IsLower ? (std::min)(cols,actual_k2) : cols - k2;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  330","line":"      // size of the triangular part"},
{"lineNum":"  331","line":"      Index ts = (IsLower && actual_k2>=cols) ? 0 : actual_kc;"},
{"lineNum":"  332","line":""},
{"lineNum":"  333","line":"      Scalar* geb = blockB+ts*ts;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  334","line":"      geb = geb + internal::first_aligned<PacketBytes>(geb,PacketBytes/sizeof(Scalar));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  335","line":""},
{"lineNum":"  336","line":"      pack_rhs(geb, rhs.getSubMapper(actual_k2,IsLower ? 0 : k2), actual_kc, rs);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  337","line":""},
{"lineNum":"  338","line":"      // pack the triangular part of the rhs padding the unrolled blocks with zeros"},
{"lineNum":"  339","line":"      if(ts>0)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  340","line":"      {"},
{"lineNum":"  341","line":"        for (Index j2=0; j2<actual_kc; j2+=SmallPanelWidth)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  342","line":"        {"},
{"lineNum":"  343","line":"          Index actualPanelWidth = std::min<Index>(actual_kc-j2, SmallPanelWidth);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  344","line":"          Index actual_j2 = actual_k2 + j2;","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  345","line":"          Index panelOffset = IsLower ? j2+actualPanelWidth : 0;"},
{"lineNum":"  346","line":"          Index panelLength = IsLower ? actual_kc-j2-actualPanelWidth : j2;"},
{"lineNum":"  347","line":"          // general part"},
{"lineNum":"  348","line":"          pack_rhs_panel(blockB+j2*actual_kc,","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  349","line":"                         rhs.getSubMapper(actual_k2+panelOffset, actual_j2),","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  350","line":"                         panelLength, actualPanelWidth,"},
{"lineNum":"  351","line":"                         actual_kc, panelOffset);"},
{"lineNum":"  352","line":""},
{"lineNum":"  353","line":"          // append the triangular part via a temporary buffer"},
{"lineNum":"  354","line":"          for (Index j=0;j<actualPanelWidth;++j)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  355","line":"          {"},
{"lineNum":"  356","line":"            if (SetDiag)"},
{"lineNum":"  357","line":"              triangularBuffer.coeffRef(j,j) = rhs(actual_j2+j,actual_j2+j);"},
{"lineNum":"  358","line":"            for (Index k=IsLower ? j+1 : 0; IsLower ? k<actualPanelWidth : k<j; ++k)"},
{"lineNum":"  359","line":"              triangularBuffer.coeffRef(k,j) = rhs(actual_j2+k,actual_j2+j);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  360","line":"          }"},
{"lineNum":"  361","line":""},
{"lineNum":"  362","line":"          pack_rhs_panel(blockB+j2*actual_kc,","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  363","line":"                         RhsMapper(triangularBuffer.data(), triangularBuffer.outerStride()),"},
{"lineNum":"  364","line":"                         actualPanelWidth, actualPanelWidth,"},
{"lineNum":"  365","line":"                         actual_kc, j2);"},
{"lineNum":"  366","line":"        }"},
{"lineNum":"  367","line":"      }"},
{"lineNum":"  368","line":""},
{"lineNum":"  369","line":"      for (Index i2=0; i2<rows; i2+=mc)","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  370","line":"      {"},
{"lineNum":"  371","line":"        const Index actual_mc = (std::min)(mc,rows-i2);","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  372","line":"        pack_lhs(blockA, lhs.getSubMapper(i2, actual_k2), actual_kc, actual_mc);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  373","line":""},
{"lineNum":"  374","line":"        // triangular kernel"},
{"lineNum":"  375","line":"        if(ts>0)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  376","line":"        {"},
{"lineNum":"  377","line":"          for (Index j2=0; j2<actual_kc; j2+=SmallPanelWidth)","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  378","line":"          {"},
{"lineNum":"  379","line":"            Index actualPanelWidth = std::min<Index>(actual_kc-j2, SmallPanelWidth);","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  380","line":"            Index panelLength = IsLower ? actual_kc-j2 : j2+actualPanelWidth;"},
{"lineNum":"  381","line":"            Index blockOffset = IsLower ? j2 : 0;"},
{"lineNum":"  382","line":""},
{"lineNum":"  383","line":"            gebp_kernel(res.getSubMapper(i2, actual_k2 + j2),","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  384","line":"                        blockA, blockB+j2*actual_kc,"},
{"lineNum":"  385","line":"                        actual_mc, panelLength, actualPanelWidth,"},
{"lineNum":"  386","line":"                        alpha,","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  387","line":"                        actual_kc, actual_kc,  // strides"},
{"lineNum":"  388","line":"                        blockOffset, blockOffset);// offsets"},
{"lineNum":"  389","line":"          }"},
{"lineNum":"  390","line":"        }"},
{"lineNum":"  391","line":"        gebp_kernel(res.getSubMapper(i2, IsLower ? 0 : k2),","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  392","line":"                    blockA, geb, actual_mc, actual_kc, rs,"},
{"lineNum":"  393","line":"                    alpha,","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  394","line":"                    -1, -1, 0, 0);"},
{"lineNum":"  395","line":"      }"},
{"lineNum":"  396","line":"    }"},
{"lineNum":"  397","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  398","line":""},
{"lineNum":"  399","line":"/***************************************************************************"},
{"lineNum":"  400","line":"* Wrapper to product_triangular_matrix_matrix"},
{"lineNum":"  401","line":"***************************************************************************/"},
{"lineNum":"  402","line":""},
{"lineNum":"  403","line":"} // end namespace internal"},
{"lineNum":"  404","line":""},
{"lineNum":"  405","line":"namespace internal {"},
{"lineNum":"  406","line":"template<int Mode, bool LhsIsTriangular, typename Lhs, typename Rhs>"},
{"lineNum":"  407","line":"struct triangular_product_impl<Mode,LhsIsTriangular,Lhs,false,Rhs,false>"},
{"lineNum":"  408","line":"{"},
{"lineNum":"  409","line":"  template<typename Dest> static void run(Dest& dst, const Lhs &a_lhs, const Rhs &a_rhs, const typename Dest::Scalar& alpha)"},
{"lineNum":"  410","line":"  {","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  411","line":"    typedef typename Lhs::Scalar  LhsScalar;"},
{"lineNum":"  412","line":"    typedef typename Rhs::Scalar  RhsScalar;"},
{"lineNum":"  413","line":"    typedef typename Dest::Scalar Scalar;"},
{"lineNum":"  414","line":""},
{"lineNum":"  415","line":"    typedef internal::blas_traits<Lhs> LhsBlasTraits;"},
{"lineNum":"  416","line":"    typedef typename LhsBlasTraits::DirectLinearAccessType ActualLhsType;"},
{"lineNum":"  417","line":"    typedef typename internal::remove_all<ActualLhsType>::type ActualLhsTypeCleaned;"},
{"lineNum":"  418","line":"    typedef internal::blas_traits<Rhs> RhsBlasTraits;"},
{"lineNum":"  419","line":"    typedef typename RhsBlasTraits::DirectLinearAccessType ActualRhsType;"},
{"lineNum":"  420","line":"    typedef typename internal::remove_all<ActualRhsType>::type ActualRhsTypeCleaned;"},
{"lineNum":"  421","line":""},
{"lineNum":"  422","line":"    typename internal::add_const_on_value_type<ActualLhsType>::type lhs = LhsBlasTraits::extract(a_lhs);"},
{"lineNum":"  423","line":"    typename internal::add_const_on_value_type<ActualRhsType>::type rhs = RhsBlasTraits::extract(a_rhs);"},
{"lineNum":"  424","line":""},
{"lineNum":"  425","line":"    LhsScalar lhs_alpha = LhsBlasTraits::extractScalarFactor(a_lhs);"},
{"lineNum":"  426","line":"    RhsScalar rhs_alpha = RhsBlasTraits::extractScalarFactor(a_rhs);"},
{"lineNum":"  427","line":"    Scalar actualAlpha = alpha * lhs_alpha * rhs_alpha;","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  428","line":""},
{"lineNum":"  429","line":"    typedef internal::gemm_blocking_space<(Dest::Flags&RowMajorBit) ? RowMajor : ColMajor,Scalar,Scalar,"},
{"lineNum":"  430","line":"              Lhs::MaxRowsAtCompileTime, Rhs::MaxColsAtCompileTime, Lhs::MaxColsAtCompileTime,4> BlockingType;"},
{"lineNum":"  431","line":""},
{"lineNum":"  432","line":"    enum { IsLower = (Mode&Lower) == Lower };"},
{"lineNum":"  433","line":"    Index stripedRows  = ((!LhsIsTriangular) || (IsLower))  ? lhs.rows() : (std::min)(lhs.rows(),lhs.cols());","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  434","line":"    Index stripedCols  = ((LhsIsTriangular)  || (!IsLower)) ? rhs.cols() : (std::min)(rhs.cols(),rhs.rows());"},
{"lineNum":"  435","line":"    Index stripedDepth = LhsIsTriangular ? ((!IsLower) ? lhs.cols() : (std::min)(lhs.cols(),lhs.rows()))","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  436","line":"                                         : ((IsLower)  ? rhs.rows() : (std::min)(rhs.rows(),rhs.cols()));"},
{"lineNum":"  437","line":""},
{"lineNum":"  438","line":"    BlockingType blocking(stripedRows, stripedCols, stripedDepth, 1, false);"},
{"lineNum":"  439","line":""},
{"lineNum":"  440","line":"    internal::product_triangular_matrix_matrix<Scalar, Index,","class":"lineNoCov","hits":"0","possible_hits":"7",},
{"lineNum":"  441","line":"      Mode, LhsIsTriangular,"},
{"lineNum":"  442","line":"      (internal::traits<ActualLhsTypeCleaned>::Flags&RowMajorBit) ? RowMajor : ColMajor, LhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  443","line":"      (internal::traits<ActualRhsTypeCleaned>::Flags&RowMajorBit) ? RowMajor : ColMajor, RhsBlasTraits::NeedToConjugate,"},
{"lineNum":"  444","line":"      (internal::traits<Dest          >::Flags&RowMajorBit) ? RowMajor : ColMajor, Dest::InnerStrideAtCompileTime>"},
{"lineNum":"  445","line":"      ::run("},
{"lineNum":"  446","line":"        stripedRows, stripedCols, stripedDepth,   // sizes"},
{"lineNum":"  447","line":"        &lhs.coeffRef(0,0), lhs.outerStride(),    // lhs info"},
{"lineNum":"  448","line":"        &rhs.coeffRef(0,0), rhs.outerStride(),    // rhs info"},
{"lineNum":"  449","line":"        &dst.coeffRef(0,0), dst.innerStride(), dst.outerStride(),    // result info"},
{"lineNum":"  450","line":"        actualAlpha, blocking"},
{"lineNum":"  451","line":"      );"},
{"lineNum":"  452","line":""},
{"lineNum":"  453","line":"    // Apply correction if the diagonal is unit and a scalar factor was nested:"},
{"lineNum":"  454","line":"    if ((Mode&UnitDiag)==UnitDiag)"},
{"lineNum":"  455","line":"    {"},
{"lineNum":"  456","line":"      if (LhsIsTriangular && lhs_alpha!=LhsScalar(1))"},
{"lineNum":"  457","line":"      {"},
{"lineNum":"  458","line":"        Index diagSize = (std::min)(lhs.rows(),lhs.cols());"},
{"lineNum":"  459","line":"        dst.topRows(diagSize) -= ((lhs_alpha-LhsScalar(1))*a_rhs).topRows(diagSize);"},
{"lineNum":"  460","line":"      }"},
{"lineNum":"  461","line":"      else if ((!LhsIsTriangular) && rhs_alpha!=RhsScalar(1))"},
{"lineNum":"  462","line":"      {"},
{"lineNum":"  463","line":"        Index diagSize = (std::min)(rhs.rows(),rhs.cols());"},
{"lineNum":"  464","line":"        dst.leftCols(diagSize) -= (rhs_alpha-RhsScalar(1))*a_lhs.leftCols(diagSize);"},
{"lineNum":"  465","line":"      }"},
{"lineNum":"  466","line":"    }"},
{"lineNum":"  467","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"8",},
{"lineNum":"  468","line":"};"},
{"lineNum":"  469","line":""},
{"lineNum":"  470","line":"} // end namespace internal"},
{"lineNum":"  471","line":""},
{"lineNum":"  472","line":"} // end namespace Eigen"},
{"lineNum":"  473","line":""},
{"lineNum":"  474","line":"#endif // EIGEN_TRIANGULAR_MATRIX_MATRIX_H"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 79, "covered" : 0,};
var merged_data = [];
