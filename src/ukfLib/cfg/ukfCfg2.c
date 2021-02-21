/**
 * @file ukfCfg2.c
 * @brief Config file for UKF example from 
 * https://thepoorengineer.com/en/ekf-impl/
 * @version 0.1
 * @date 2021-02-20
 * 
 * State vector
 * x[0] = q0
 * x[1] = q1
 * x[2] = q2
 * x[2] = q3
 * x[3] = bg1
 * x[4] = bg2
 * x[5] = bg3
 * 
 * 
 * State transition matrix
 * {1,    0,    0,    0,  T/2*q1,  T/2*q2,  T/2*q3}
 * {0,    1,    0,    0, -T/2*q0,  T/2*q3, -T/2*q2}
 * {0,    0,    1,    0, -T/2*q3, -T/2*q0,  T/2*q1}
 * {0,    0,    0,    1,  T/2*q2, -T/2*q1, -T/2*q0}
 * {0,    0,    0,    0,       1,       0,       0}
 * {0,    0,    0,    0,       0,       1,       0}
 * {0,    0,    0,    0,       0,       0,       1}
 * 
 */

#include "ukfCfg2.h"
#include <stdint.h>
#include <math.h>

static void Fx0(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx1(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx2(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx3(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx4(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx5(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx6(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);

static void Hy0(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx);
static void Hy1(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx);

static tPredictFcn PredictFcn[Lx] = {&Fx0, &Fx1, &Fx2, &Fx3, &Fx4, &Fx5, &Fx6};
static tObservFcn   ObservFcn[Ly] = {&Hy0, &Hy1};

//! UKF Processing matrix
static float Sc_vector[1][3] = {{1, 2, 0}};
static float Wm_weight_vector[1][2*Lx + 1] = {0};
static float Wc_weight_vector[1][2*Lx + 1] = {0};
//static float u_system_input[4][1] = {{0},{0},{0},{0}};
//static float u_prev_system_input[4][1] = {{0},{0},{0},{0}};
static float y_meas[Ly][1] = {0};
static float y_predicted_mean[Ly][1] = {0};
static float x_system_states[Lx][1] = {{0}, {0}, {50}, {50}};
static float x_system_states_ic[Lx][1] = {{0}, {0}, {50}, {50}};
static float x_system_states_correction[Lx][1] = {0};
static float X_sigma_points[Lx][2*Lx + 1] = {0};

//! Sigma points Y(k|k-1) = y_m
static float Y_sigma_points[Ly][2*Lx + 1] = {0};

//! State covariance  P(k|k-1) = P_m, P(k)= P
static float Pxx_error_covariance[Lx][Lx] = {0};

static float Pxx0_init_error_covariance[Lx][Lx] = {
        {1, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 1},
};

//! Process noise covariance Q : initial noise assumptions
static float Qxx_process_noise_cov[Lx][Lx] = {
        {1, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 1},
};

//! Output noise covariance: initial noise assumptions
static float Ryy0_init_out_covariance[Ly][Ly] =
    {
        /*  y1, y2         */
        {1, 0}, /* y1 */
        {0, 1}, /* y2 */
};

//! Output covariance Pyy = R (initial assumption)
static float Pyy_out_covariance[Ly][Ly] = {0};

static float Pyy_out_covariance_copy[Ly][Ly] = {0};

//! cross-covariance of state and output
static float Pxy_cross_covariance[Lx][Ly] = {0};

//! Kalman gain matrix
static float K_kalman_gain[Lx][Ly] = {0};

static float Pxx_covariance_correction[Lx][Lx] = {0};

static float I_identity_matrix[Ly][Ly] = {0};

tUkfMatrix UkfMatrixCfg2 = {
    .Sc_vector                      = {NROWS(Sc_vector),        NCOL(Sc_vector),        &Sc_vector[0][0]},
    .Wm_weight_vector               = {NROWS(Wm_weight_vector), NCOL(Wm_weight_vector), &Wm_weight_vector[0][0]},
    .Wc_weight_vector               = {NROWS(Wc_weight_vector), NCOL(Wc_weight_vector), &Wc_weight_vector[0][0]},
    .x_system_states                = {NROWS(x_system_states),  NCOL(x_system_states),  &x_system_states[0][0]},
    .x_system_states_ic             = {NROWS(x_system_states_ic), NCOL(x_system_states_ic), &x_system_states_ic[0][0]},
    .x_system_states_limits         = {0, 0, NULL},
    .x_system_states_limits_enable  = {0, 0, NULL},
    .x_system_states_correction     = {NROWS(x_system_states_correction), NCOL(x_system_states_correction), &x_system_states_correction[0][0]},
    .u_system_input                 = {0, 0, NULL},
    .u_prev_system_input            = {0, 0, NULL},
    .X_sigma_points                 = {NROWS(X_sigma_points), NCOL(X_sigma_points), &X_sigma_points[0][0]},
    .Y_sigma_points                 = {NROWS(Y_sigma_points), NCOL(Y_sigma_points), &Y_sigma_points[0][0]},
    .y_predicted_mean               = {NROWS(y_predicted_mean), NCOL(y_predicted_mean), &y_predicted_mean[0][0]},
    .y_meas                         = {NROWS(y_meas), NCOL(y_meas), &y_meas[0][0]},
    .Pyy_out_covariance             = {NROWS(Pyy_out_covariance), NCOL(Pyy_out_covariance), &Pyy_out_covariance[0][0]},
    .Pyy_out_covariance_copy        = {NROWS(Pyy_out_covariance_copy), NCOL(Pyy_out_covariance_copy), &Pyy_out_covariance_copy[0][0]},
    .Ryy0_init_out_covariance       = {NROWS(Ryy0_init_out_covariance), NCOL(Ryy0_init_out_covariance), &Ryy0_init_out_covariance[0][0]},
    .Pxy_cross_covariance           = {NROWS(Pxy_cross_covariance), NCOL(Pxy_cross_covariance), &Pxy_cross_covariance[0][0]},
    .Pxx_error_covariance           = {NROWS(Pxx_error_covariance), NCOL(Pxx_error_covariance), &Pxx_error_covariance[0][0]},
    .Pxx0_init_error_covariance     = {NROWS(Pxx0_init_error_covariance), NCOL(Pxx0_init_error_covariance), &Pxx0_init_error_covariance[0][0]},
    .Qxx_process_noise_cov          = {NROWS(Qxx_process_noise_cov), NCOL(Qxx_process_noise_cov), &Qxx_process_noise_cov[0][0]},
    .K_kalman_gain                  = {NROWS(K_kalman_gain), NCOL(K_kalman_gain), &K_kalman_gain[0][0]},
    .I_identity_matrix              = {NROWS(I_identity_matrix), NCOL(I_identity_matrix), &I_identity_matrix[0][0]},
    .Pxx_covariance_correction      = {NROWS(Pxx_covariance_correction), NCOL(Pxx_covariance_correction), &Pxx_covariance_correction[0][0]},
    .fcnPredict                     = &PredictFcn[0],
    .fcnObserve                     = &ObservFcn[0],
    .dT                             = 0.1F
};

void Fx0(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;

    X_m->val[nCol * 0 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
}

void Fx1(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 1 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
}

void Fx2(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 2 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
    dT = dT;
}

void Fx3(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 3 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
    dT = dT;
}

void Fx4(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 4 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
    dT = dT;
}

void Fx5(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 5 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
    dT = dT;
}

void Fx6(tMatrix* u_p, tMatrix* X_p, tMatrix* X_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = X_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    X_m->val[nCol * 6 + sigmaIdx] = 
        X_p->val[nCol * 0 + sigmaIdx] + 
        X_p->val[nCol * 1 + sigmaIdx] + 
        X_p->val[nCol * 2 + sigmaIdx] + 
        X_p->val[nCol * 3 + sigmaIdx] + 
        X_p->val[nCol * 4 + sigmaIdx] + 
        X_p->val[nCol * 5 + sigmaIdx] + 
        X_p->val[nCol * 6 + sigmaIdx];

    u_p = u_p;
    dT = dT;
}

void Hy0(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx) {
    static const float N1 = 20;
    static const float E1 = 0;
    float term1;
    float term2;
    const uint8_t nCol = pY_m->ncol;

    term1 = pX_m->val[nCol * 0 + sigmaIdx] - N1;
    term1 *= term1;

    term2 = pX_m->val[nCol * 1 + sigmaIdx] - E1;
    term2 *= term2;

    pY_m->val[sigmaIdx] = sqrtf(term1 + term2);

    pu = pu;
}

void Hy1(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx) {
    static const float N2 = 0;
    static const float E2 = 20;
    float term1;
    float term2;
    const uint8_t nCol = pY_m->ncol;

    term1 = pX_m->val[nCol * 0 + sigmaIdx] - N2;
    term1 *= term1;

    term2 = pX_m->val[nCol * 1 + sigmaIdx] - E2;
    term2 *= term2;

    pY_m->val[nCol * 1 + sigmaIdx] = sqrtf(term1 + term2);

    pu = pu;
}
