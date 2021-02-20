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
    uint8_t nrow;
    uint8_t ncol;
    float* val;
} tMatrix;

typedef struct sMatrixBool {
    uint8_t nrow;
    uint8_t ncol;
    uint8_t* val;
} tMatrixBool;

mtxResultInfo mtx_diagsum   (tMatrix *pSrc, float *diagsum);
mtxResultInfo mtx_transp_square (tMatrix *pSrc);
mtxResultInfo mtx_transp_dest   (const tMatrix *pSrc, tMatrix *pDst);
mtxResultInfo mtx_mul           (const tMatrix *pSrc1, const tMatrix *pSrc2, tMatrix *pDst);
mtxResultInfo mtx_mul_src2tr    (const tMatrix *pSrc1, const tMatrix *pSrc2, tMatrix *pDst);
mtxResultInfo mtx_chol_lower    (tMatrix *pSrc);
mtxResultInfo mtx_chol_upper    (tMatrix *pSrc);
mtxResultInfo mtx_inv   (tMatrix *pSrc, tMatrix *pDst);
mtxResultInfo mtx_add   (tMatrix *pDst, const tMatrix *pSrc);
mtxResultInfo mtx_sub   (tMatrix *pDst, const tMatrix *pSrc);
mtxResultInfo mtx_mul_scalar    (tMatrix *pSrc, float scalar);
mtxResultInfo mtx_sub_scalar    (tMatrix *pSrc, float scalar);
mtxResultInfo mtx_add_scalar    (tMatrix *pSrc, float scalar);
mtxResultInfo mtx_cpy           (tMatrix *pDst, const tMatrix *pSrc);
mtxResultInfo mtx_identity      (tMatrix *pSrc);
mtxResultInfo mtx_zeros         (tMatrix *pSrc);
mtxResultInfo mtx_print         (const tMatrix *A);

#endif
