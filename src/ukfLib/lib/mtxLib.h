/**
 * @file mtxLib.h
 * @brief Basic matrix operations header.
 * @version 0.1
 * @date 2021-02-20
 */

#ifndef MTXLIB_FILE
#define MTXLIB_FILE

#include <math.h>
#include <stdint.h>

//! Macros definiton
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
    uint16_t nelem;
    uint8_t nrow;
    uint8_t ncol;
    float* val;
} tMatrix;

typedef struct sMatrixBool {
    uint16_t nelem;
    uint8_t nrow;
    uint8_t ncol;
    uint8_t* val;
} tMatrixBool;

mtxResultInfo mtx_init_bool(tMatrixBool* const pSrc, uint8_t* const pValue, const uint8_t nrow, const uint8_t ncol, const uint16_t nelem);
mtxResultInfo mtx_init(tMatrix* const pSrc, float* const pValue, const uint8_t nrow, const uint8_t ncol, const uint16_t nelem);
mtxResultInfo mtx_mul(tMatrix const* const pSrc1, tMatrix const* const pSrc2, tMatrix* const pDst);
mtxResultInfo mtx_transp_square(tMatrix* const pSrc);
mtxResultInfo mtx_transp_dest(tMatrix const* const pSrc, tMatrix* const pDst);
mtxResultInfo mtx_diagsum(tMatrix* pSrc, float* diagsum);
mtxResultInfo mtx_chol_upper(tMatrix* const pSrc);
mtxResultInfo mtx_chol_lower(tMatrix* const pSrc);
mtxResultInfo mtx_inv(tMatrix* const pSrc, tMatrix* const pDst);
mtxResultInfo mtx_add(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_sub(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_mul_scalar(tMatrix* const pSrc, const float scalar);
mtxResultInfo mtx_add_scalar(tMatrix* const pSrc, const float scalar);
mtxResultInfo mtx_sub_scalar(tMatrix* const pSrc, const float scalar);
mtxResultInfo mtx_cpy(tMatrix* const pDst, tMatrix const* const pSrc);
mtxResultInfo mtx_identity(tMatrix* const pSrc);
mtxResultInfo mtx_zeros(tMatrix* const pSrc);
mtxResultInfo mtx_mul_src2tr(tMatrix const* const pSrc1, tMatrix const* const pSrc2, tMatrix* const pDst);
mtxResultInfo mtx_print(tMatrix const *A);

#endif
