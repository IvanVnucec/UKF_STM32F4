/**
 * @file mtxLib.c
 * @brief Matrix operations
 * @version 0.1
 * @date 2021-02-20
 */

#include <stdio.h>
#include <stdint.h>
#include "mtxLib.h"

/**
 * @brief 
 * 
 * @param pSrc 
 * @param diagsum 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_diagsum(tMatrix *pSrc, float *diagsum) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float const *const pSrcL = (float *)pSrc->val;
    const uint8_t ncol = pSrc->ncol;
    uint16_t eIdx;
    float sum = pSrcL[0];
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    if (pSrc->nrow == ncol) {
        for (eIdx = 1; eIdx < nelem; eIdx++) {
            const uint16_t cmpLeft = (uint16_t)(eIdx / ncol);

            sum += eIdx < ncol ? 0 : cmpLeft == eIdx % (cmpLeft * ncol) ? pSrcL[eIdx]
                                                                        : 0;
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    *diagsum = sum;

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_transp_square(tMatrix *pSrc) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    const uint8_t nrow = pSrc->nrow;
    const uint8_t ncol = pSrc->ncol;
    float *const pSrcL = (float *)pSrc->val;
    uint8_t row, col;
    float temp;

    if (nrow == ncol) {
        for (row = 0; row < nrow; row++) {
            for (col = 0; col < ncol; col++) {
                if (row != col && row < col) {
                    temp = pSrcL[nrow * row + col];
                    pSrcL[ncol * row + col] = pSrcL[ncol * col + row];
                    pSrcL[ncol * col + row] = temp;
                }
            }
        }
    } else {
        ResultL = MTX_NOT_SQUARE;
    }

    return ResultL;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @param pDst 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_transp_dest(const tMatrix *pSrc, tMatrix *pDst) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float const *const pSrcL = (float *)pSrc->val;
    float *const pDstL = (float *)pDst->val;
    const uint8_t nRowSrcL = pSrc->nrow;
    const uint8_t nColSrcL = pSrc->ncol;
    const uint8_t nRowDstL = pDst->nrow;
    const uint8_t nColDstL = pDst->ncol;
    uint8_t row, col;

    if (nRowSrcL == nColDstL || nColSrcL == nRowDstL) {
        for (row = 0; row < nRowDstL; row++) {
            for (col = 0; col < nColDstL; col++) {
                pDstL[nColDstL * row + col] = pSrcL[nColSrcL * col + row];
            }
        }
    } else {
        ResultL = MTX_SIZE_MISMATCH;
    }

    return ResultL;
}

/**
 * @brief 
 * 
 * @param pSrc1 
 * @param pSrc2 
 * @param pDst 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_mul(const tMatrix *pSrc1, const tMatrix *pSrc2, tMatrix *pDst) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float const *const pSrc1L = (float *)pSrc1->val;
    float const *const pSrc2L = (float *)pSrc2->val;
    float *const pDstL = (float *)pDst->val;
    uint8_t row, col, k;
    float sum;

    if (pSrc1->ncol == pSrc2->nrow) {
        for (row = 0; row < pSrc1->nrow; row++) {
            for (col = 0; col < pSrc2->ncol; col++) {
                sum = 0;
                for (k = 0; k < pSrc1->ncol; k++) {
                    sum += pSrc1L[pSrc1->ncol * row + k] * pSrc2L[pSrc2->ncol * k + col];
                }
                pDstL[pDst->ncol * row + col] = sum;
            }
        }
    } else {
        ResultL = MTX_SIZE_MISMATCH;
    }

    return ResultL;
}

/**
 * @brief Special multiplication Dst=Src1*Src2'
 * Special function for multiplication of Src1 matrix with 
 * transpose image of Src2 matrix. Be sure that array size are
 * suitable for multiplication after Src2 transpose  
 * 
 * @param pSrc1 
 * @param pSrc2 
 * @param pDst 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_mul_src2tr(const tMatrix *pSrc1, const tMatrix *pSrc2, tMatrix *pDst) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float const *const pSrc1L = (float *)pSrc1->val;
    float const *const pSrc2L = (float *)pSrc2->val;
    float *const pDstL = (float *)pDst->val;
    uint8_t rowSrc1, rowSrc2, k;
    float sum;

    if (pSrc1->ncol == pSrc2->ncol) {
        for (rowSrc1 = 0; rowSrc1 < pSrc1->nrow; rowSrc1++) {
            for (rowSrc2 = 0; rowSrc2 < pSrc2->nrow; rowSrc2++) {
                sum = 0;
                for (k = 0; k < pSrc1->ncol; k++) {
                    sum += pSrc1L[pSrc1->ncol * rowSrc1 + k] * pSrc2L[pSrc2->ncol * rowSrc2 + k];
                }
                pDstL[pDst->ncol * rowSrc1 + rowSrc2] = sum;
            }
        }
    } else {
        ResultL = MTX_SIZE_MISMATCH;
    }

    return ResultL;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_chol_lower(tMatrix *pSrc) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float *const pSrcL = pSrc->val;
    const uint8_t nrow = pSrc->nrow;
    const uint8_t ncol = pSrc->ncol;
    uint8_t col, row;
    int8_t tmp;
    float sum = 0;

    if (ncol == nrow) {
        const uint8_t mtxSize = nrow;

        for (col = 0; col < mtxSize; col++) {
            for (row = 0; row < mtxSize; row++) {
                sum = pSrcL[mtxSize * col + row];

                for (tmp = (int8_t)(col - 1); tmp >= 0; tmp--) {
                    sum -= pSrcL[mtxSize * row + tmp] * pSrcL[mtxSize * col + tmp];
                }

                pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row > col) ? (sum / pSrcL[ncol * col + col])
                                                                                 : 0;

                if ((row == col) && (sum <= 0)) {
                    ResultL = MTX_NOT_POS_DEFINED;
                }
            }
        }
    } else {
        ResultL = MTX_NOT_SQUARE;
    }

    return ResultL;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_chol_upper(tMatrix *pSrc) {
    mtxResultInfo ResultL = MTX_OPERATION_OK;
    float *const pSrcL = pSrc->val;
    const uint8_t nrow = pSrc->nrow;
    const uint8_t ncol = pSrc->ncol;
    uint8_t col, row;
    int8_t tmp;
    float sum = 0;

    if (ncol == nrow) {
        for (row = 0; row < nrow; row++) {
            for (col = 0; col < ncol; col++) {
                sum = pSrcL[ncol * row + col];

                for (tmp = (int8_t)(row - 1); tmp >= 0; tmp--)  // tmp could be calc negative
                {
                    sum -= pSrcL[ncol * tmp + row] * pSrcL[ncol * tmp + col];
                }

                pSrcL[ncol * row + col] = (row == col) ? sqrt(sum) : (row < col) ? (sum / pSrcL[ncol * row + row])
                                                                                 : 0;

                if ((row == col) && (sum <= 0)) {
                    ResultL = MTX_NOT_POS_DEFINED;
                }
            }
        }
    } else {
        ResultL = MTX_NOT_SQUARE;
    }

    return ResultL;
}

/**
 * @brief Matrix inverse.
 * 
 * @param pSrc Square matrix
 * @param pDst At the begining should point to identity matrix!!
 * @note  
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_inv(tMatrix *pSrc, tMatrix *pDst) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    const uint8_t nrow = pSrc->nrow;
    const uint8_t ncol = pSrc->ncol;
    uint8_t j, i;
    uint8_t k = 0;
    uint8_t l = 0;
    float s = 0;
    float t = 0;

    /* TODO: implement pDst so it doesnt need to be initialized to an identity 
    matrix beforehand, rather in this function IvanVnucec 
    */

    if (nrow == ncol) {
        for (j = 0; j < nrow; j++) {
            for (i = j; i < nrow; i++) {
                if (0 != pSrc->val[ncol * i + j]) {
                    for (k = 0; k < nrow; k++) {
                        s = pSrc->val[ncol * j + k];
                        pSrc->val[ncol * j + k] = pSrc->val[ncol * i + k];
                        pSrc->val[ncol * i + k] = s;

                        s = pDst->val[ncol * j + k];
                        pDst->val[ncol * j + k] = pDst->val[ncol * i + k];
                        pDst->val[ncol * i + k] = s;
                    }

                    t = 1 / pSrc->val[ncol * j + j];

                    for (k = 0; k < nrow; k++) {
                        pSrc->val[ncol * j + k] = t * pSrc->val[ncol * j + k];
                        pDst->val[ncol * j + k] = t * pDst->val[ncol * j + k];
                    }

                    for (l = 0; l < nrow; l++) {
                        if (l != j) {
                            t = -pSrc->val[ncol * l + j];
                            for (k = 0; k < nrow; k++) {
                                pSrc->val[ncol * l + k] += t * pSrc->val[ncol * j + k];
                                pDst->val[ncol * l + k] += t * pDst->val[ncol * j + k];
                            }
                        }
                    }
                }
                break;
            }
            if (0 == pSrc->val[ncol * l + k]) {
                Result = MTX_SINGULAR;
            }
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pDst 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_add(tMatrix *pDst, const tMatrix *pSrc) {
    uint8_t Result = MTX_OPERATION_OK;
    float *const pDstL = (float *)pDst->val;
    float const *const pSrcL = (float *)pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow) {
        for (eIdx = 0; eIdx < nelem; eIdx++) {
            pDstL[eIdx] += pSrcL[eIdx];
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pDst 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_sub(tMatrix *pDst, const tMatrix *pSrc) {
    uint8_t Result = MTX_OPERATION_OK;
    float *const pDstL = (float *)pDst->val;
    float const *const pSrcL = (float *)pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow) {
        for (eIdx = 0; eIdx < nelem; eIdx++) {
            pDstL[eIdx] -= pSrcL[eIdx];
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @param scalar 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_mul_scalar(tMatrix *pSrc, float scalar) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDst = pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    for (eIdx = 0; eIdx < nelem; eIdx++) {
        pDst[eIdx] *= scalar;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @param scalar 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_sub_scalar(tMatrix *pSrc, float scalar) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDst = pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    for (eIdx = 0; eIdx < nelem; eIdx++) {
        pDst[eIdx] -= scalar;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @param scalar 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_add_scalar(tMatrix *pSrc, float scalar) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDst = pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    for (eIdx = 0; eIdx < nelem; eIdx++) {
        pDst[eIdx] += scalar;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pDst 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_cpy(tMatrix *pDst, const tMatrix *pSrc) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDstL = pDst->val;
    float const *const pSrcL = pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    if (pDst->ncol == pSrc->ncol && pDst->nrow == pSrc->nrow) {
        for (eIdx = 0; eIdx < nelem; eIdx++) {
            pDstL[eIdx] = pSrcL[eIdx];
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_identity(tMatrix *pSrc) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDst = (float *)pSrc->val;
    const uint8_t nCol = pSrc->ncol;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    if (pSrc->nrow == nCol) {
        pDst[0] = 1;

        for (eIdx = 1; eIdx < nelem; eIdx++) {
            const uint16_t cmpLeft = (uint16_t)(eIdx / nCol);

            /* TODO: Optimize this so we initialize matrix to all zeros and then with
             * the for loop only initialize diagonals to 1.0 */
            pDst[eIdx] = eIdx < nCol ? 0 : cmpLeft == eIdx % (cmpLeft * nCol) ? 1
                                                                              : 0;
        }
    } else {
        Result = MTX_SIZE_MISMATCH;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param pSrc 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_zeros(tMatrix *pSrc) {
    mtxResultInfo Result = MTX_OPERATION_OK;
    float *const pDst = (float *)pSrc->val;
    uint16_t eIdx;
    const uint16_t nelem = pSrc->ncol * pSrc->nrow;

    for (eIdx = 0; eIdx < nelem; eIdx++) {
        pDst[eIdx] = 0;
    }

    return Result;
}

/**
 * @brief 
 * 
 * @param A 
 * @return mtxResultInfo 
 */
mtxResultInfo mtx_print(const tMatrix *A) {
    int i, j;
    mtxResultInfo Result = MTX_OPERATION_OK;

    for (i = 0; i < A->nrow; i++) {
        for (j = 0; j < A->ncol; j++) {
            printf("% -3.5f ", A->val[A->ncol * i + j]);
        }
        printf("\n");
    }
    printf("\n");

    return Result;
}
