/******************************************************************************************************************************************************************************************************\
 *** 
 *** Description       : IMPLEMENTATION OF BASIC MATRIX OPERATION
 *** Codefile          : mtxLib.h
 ***
 *** MIT License
 ***
 *** Copyright (c) 2017 ivo-georgiev
 ***  
 *** Permission is hereby granted, free of charge, to any person obtaining a copy
 *** of this software and associated documentation files (the "Software"), to deal
 *** in the Software without restriction, including without limitation the rights
 *** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *** copies of the Software, and to permit persons to whom the Software is
 *** furnished to do so, subject to the following conditions:
 ***    
 *** The above copyright notice and this permission notice shall be included in all
 *** copies or substantial portions of the Software.
 ***      
 *** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *** OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *** SOFTWARE.
\******************************************************************************************************************************************************************************************************/
#ifndef MTXLIB_FILE
#define MTXLIB_FILE

#include <math.h>

#include "System_Types.h"
/*---------------------------------------------*/
/*         Macros definiton                    */
/*---------------------------------------------*/
#define NCOL(arr) (sizeof(arr[0]) / sizeof(arr[0][0]))
#define NROWS(arr) (sizeof(arr) / sizeof(arr[0]))
#define COLXROW(arr) (sizeof(arr) / sizeof(arr[0][0]))

#define MTX_OPERATION_OK (0UL)
#define MTX_SINGULAR (251UL)
#define MTX_SIZE_MISMATCH (252UL)
#define MTX_NOT_SQUARE (253UL)
#define MTX_NOT_POS_DEFINED (254UL)
#define MTX_OPERATION_ERROR (255UL)

typedef int mtxResultInfo;

typedef struct sMatrix {
    uint16 nelem;
    uint8 nrow;
    uint8 ncol;
    double* val;
} tMatrix;

typedef struct sMatrixBool {
    uint16 nelem;
    uint8 nrow;
    uint8 ncol;
    boolean* val;
} tMatrixBool;

mtxResultInfo mtx_init_bool(tMatrixBool* const pSrc, boolean* const pValue, const uint8 nrow, const uint8 ncol, const uint16 nelem);
mtxResultInfo mtx_init_f64(tMatrix* const pSrc, float64* const pValue, const uint8 nrow, const uint8 ncol, const uint16 nelem);
mtxResultInfo mtx_mul_f64(tMatrix const* const pSrc1, tMatrix const* const pSrc2, tMatrix* const pDst);
mtxResultInfo mtx_transp_square_f64(tMatrix* const pSrc);
mtxResultInfo mtx_transp_dest_f64(tMatrix const* const pSrc, tMatrix* const pDst);
mtxResultInfo mtx_diagsum_f64(tMatrix* pSrc, double* diagsum);
mtxResultInfo mtx_chol_upper_f64(tMatrix* const pSrc);
mtxResultInfo mtx_chol_lower_f64(tMatrix* const pSrc);
mtxResultInfo mtx_inv_f64(tMatrix* const pSrc, tMatrix* const pDst);
mtxResultInfo mtx_add_f64(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_sub_f64(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_mul_scalar_f64(tMatrix* const pSrc, const float64 scalar);
mtxResultInfo mtx_add_scalar_f64(tMatrix* const pSrc, const float64 scalar);
mtxResultInfo mtx_sub_scalar_f64(tMatrix* const pSrc, const float64 scalar);
mtxResultInfo mtx_cpy_f64(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_identity_f64(tMatrix* const pSrc);
mtxResultInfo mtx_zeros_f64(tMatrix* const pSrc);
mtxResultInfo mtx_mul_src2tr_f64(tMatrix const* const pSrc1, tMatrix const* const pSrc2, tMatrix* const pDst);

#endif
