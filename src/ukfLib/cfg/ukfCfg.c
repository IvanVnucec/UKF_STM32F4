/**
 * @file ukfCfg.c
 * @brief Config file for UKF example from 
 * https://yugu.faculty.wvu.edu/files/d/2cbb566f-9936-4033-bb1c-6d887c30d45a/irl_wvu_online_ukf_implementation_v1-0_06_28_2013.pdf
 * @version 0.1
 * @date 2021-02-20
 * 
 * State vector
 * x[0] = n(k)
 * x[1] = e(k)
 * x[2] = ndot(k)
 * x[2] = edot(k)
 * 
 * 
 * State transition matrix
 * {1,   0, dT0,   0},
 * {0,   1,   0, dT0},
 * {0,   0,   1,   0},
 * {0,   0,   0,   1}
 * 
 */

#include "ukfCfg.h"
#include <stdint.h>
#include <math.h>

static void Fx1(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx2(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx3(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);
static void Fx4(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT);

static void Hy1(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx);
static void Hy2(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx);

static tPredictFcn PredictFcn[Lx] = {&Fx1, &Fx2, &Fx3, &Fx4};
static tObservFcn ObservFcn[Ly] = {&Hy1, &Hy2};

//! UKF Processing matrix
static float Sc_vector[1][3] = {{1, 2, 0}};
static float Wm_weight_vector[1][2*Lx + 1] = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
static float Wc_weight_vector[1][2*Lx + 1] = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
//static float u_system_input[4][1] = {{0},{0},{0},{0}};
//static float u_prev_system_input[4][1] = {{0},{0},{0},{0}};
static float y_meas[Ly][1] = {{0}, {0}};
static float y_predicted_mean[Ly][1] = {{0}, {0}};
static float x_system_states[Lx][1] = {{0}, {0}, {50}, {50}};
static float x_system_states_ic[Lx][1] = {{0}, {0}, {50}, {50}};
static float x_system_states_correction[Lx][1] = {{0}, {0}, {0}, {0}};
static float X_sigma_points[Lx][2*Lx + 1] =
    {
        /*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* x1 */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* x2 */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* x3 */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* x4 */
};

//! Sigma points Y(k|k-1) = y_m
static float Y_sigma_points[Ly][2*Lx + 1] =
    {
        /*  s1  s2  s3  s4  s5  s6  s7  s8  s9        */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* y1 */
        {0, 0, 0, 0, 0, 0, 0, 0, 0}, /* y2 */
};

//! State covariance  P(k|k-1) = P_m, P(k)= P
static float Pxx_error_covariance[Lx][Lx] =
    {
        /*  x1, x2, x3, x4        */
        {0, 0, 0, 0}, /* x1 */
        {0, 0, 0, 0}, /* x2 */
        {0, 0, 0, 0}, /* x3 */
        {0, 0, 0, 0}, /* x4 */
};

//State covariance initial values
/*Matthew B Rhudy : initial error covariance matrix should be defined 
based on your initialization error.
I.e., if you think your initial state is not very close, the P0 value 
should be large,
whereas if the initialization is very good (high confidence that your 
states are close to the correct values) you can assume a smaller P0 value.
I would not recommend setting P0 to zero, as this assumes there is no 
initialization error and your initial states are perfect.  
This is almost never the case.  Many times the P0 matrix is diagonal, 
with the diagonal components corresponding to the expected variance
in the corresponding state, i.e. how much deviation you might expect 
in the initialization of that state.  If you have no idea where to start, 
I recommend using an identity matrix rather than the zero matrix. */
static float Pxx0_init_error_covariance[Lx][Lx] =
    {
        /*  x1, x2, x3, x4        */
        {1, 0, 0, 0}, /* x1 */
        {0, 1, 0, 0}, /* x2 */
        {0, 0, 1, 0}, /* x3 */
        {0, 0, 0, 1}, /* x4 */
};

//! Process noise covariance Q : initial noise assumptions
/* Matthew B Rhudy : Q matrix corresponds to the uncertainty that you 
expect in your state equations.
This could include modeling errors or other uncertainties in the 
equations themselves.
Some formulations consider input measurements in the state equations 
which introduces process noise.
If you are very confident in your equations, you could set Q to zero. 
If you do that the filter will use
the noise free model to predict the state vector and will ignore any 
measurement data since your model is assumed perfect. */
static float Qxx_process_noise_cov[Lx][Lx] =
    {
        /*  x1, x2, x3, x4        */
        {0, 0, 0, 0}, /* x1 */
        {0, 0, 0, 0}, /* x2 */
        {0, 0, 4, 0}, /* x3 */
        {0, 0, 0, 4}, /* x4 */
};

//! Output noise covariance: initial noise assumptions
static float Ryy0_init_out_covariance[Ly][Ly] =
    {
        /*  y1, y2         */
        {1, 0}, /* y1 */
        {0, 1}, /* y2 */
};

//! Output covariance Pyy = R (initial assumption)
static float Pyy_out_covariance[Ly][Ly] =
    {
        /*  y1, y2         */
        {0, 0}, /* y1 */
        {0, 0}, /* y2 */
};

static float Pyy_out_covariance_copy[Ly][Ly] =
    {
        /*  y1, y2         */
        {0, 0}, /* y1 */
        {0, 0}, /* y2 */
};

//! cross-covariance of state and output
static float Pxy_cross_covariance[Lx][Ly] =
    {
        /*  y1, y2         */
        {0, 0}, /* x1 */
        {0, 0}, /* x2 */
        {0, 0}, /* x3 */
        {0, 0}, /* x4 */
};

//! Kalman gain matrix
static float K_kalman_gain[Lx][Ly] =
    {
        {0, 0},
        {0, 0},
        {0, 0},
        {0, 0},
};

static float Pxx_covariance_correction[Lx][Lx] =
    {
        /*  x1, x2, x3, x4        */
        {0, 0, 0, 0}, /* x1 */
        {0, 0, 0, 0}, /* x2 */
        {0, 0, 0, 0}, /* x3 */
        {0, 0, 0, 0}, /* x4 */
};

static float I_identity_matrix[Ly][Ly] =
    {
        {0, 0},
        {0, 0},
};

tUkfMatrix UkfMatrixCfg = {
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

/**
 * @brief Calculate predicted state 0 for each sigma point. 
 * Note  that  this  problem  has  a  linear  prediction stage 
 * X_m[0][sigmaIdx] = f(X_p, u_p) = n(k)  = n(k-1) + dT * ndot(k-1)
 * 
 * @param pu_p NULL for this system, be sure that is not used in calc
 * @param pX_p Pointer to the sigma points array at (k-1) moment 
 * @param pX_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 * @param dT Sampling time.
 */
void Fx1(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = pX_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    pX_m->val[nCol * 0 + sigmaIdx] = pX_p->val[nCol * 0 + sigmaIdx] + dT * pX_p->val[nCol * 2 + sigmaIdx];

    pu_p = pu_p;
}

/**
 * @brief Calculate predicted state 1 for each sigma point. 
 * Note  that  this  problem  has  a  linear  prediction stage 
 * X_m[1][sigmaIdx] = f(X_p, u_p) = e(k)  = e(k-1) + dT * edot(k-1)
 * 
 * @param pu_p NULL for this system, be sure that is not used in calc
 * @param pX_p Pointer to the sigma points array at (k-1) moment 
 * @param pX_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 * @param dT Sampling time.
 */
void Fx2(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = pX_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    pX_m->val[nCol * 1 + sigmaIdx] = pX_p->val[nCol * 1 + sigmaIdx] + dT * pX_p->val[nCol * 3 + sigmaIdx];

    pu_p = pu_p;
}

/**
 * @brief Calculate predicted state 2 for each sigma point. 
 * Note  that  this  problem  has  a  linear  prediction stage 
 * X_m[2][sigmaIdx] = f(X_p, u_p) = ndot(k)  = ndot(k-1)
 * 
 * @param pu_p NULL for this system, be sure that is not used in calc
 * @param pX_p Pointer to the sigma points array at (k-1) moment 
 * @param pX_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 * @param dT Sampling time.
 */
void Fx3(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = pX_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    pX_m->val[nCol * 2 + sigmaIdx] = pX_p->val[nCol * 2 + sigmaIdx];

    pu_p = pu_p;
    dT = dT;
}

/**
 * @brief  Calculate predicted state 3 for each sigma point. 
 * Note  that  this  problem  has  a  linear  prediction stage 
 * X_m[3][sigmaIdx] = f(X_p, u_p) = edot(k)  = edot(k-1)   
 * 
 * @param pu_p NULL for this system, be sure that is not used in calc
 * @param pX_p Pointer to the sigma points array at (k-1) moment 
 * @param pX_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 * @param dT Sampling time.
 */
void Fx4(tMatrix* pu_p, tMatrix* pX_p, tMatrix* pX_m, uint8_t sigmaIdx, float dT) {
    const uint8_t nCol = pX_m->ncol;  //pX_m->ncol == pX_p->ncol == 9

    pX_m->val[nCol * 3 + sigmaIdx] = pX_p->val[nCol * 3 + sigmaIdx];

    pu_p = pu_p;
    dT = dT;
}

/**
 * @brief  Calculate predicted state 3 for each sigma point. 
 * This problem has a nonlinear observation 
 * Y_m[0][sigmaIdx] = h1(X_m, u) = y1(k)  = sqrtf((n(k)-N1)^2+(e(k)-E1)^2)   
 * 
 * @param pu NULL for this system, be sure that is not used in calc
 * @param pX_m Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param pY_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 */
void Hy1(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx) {
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

/**
 * @brief Calculate predicted state 3 for each sigma point. 
 * This problem has a nonlinear observation 
 * Y_m[1][sigmaIdx] = h2(X_m[0&1][sigmaIdx], u) = y2(k)  = sqrtf( (n(k) - N2)^2 + (e(k) - E2)^2 )    
 * 
 * @param pu NULL for this system, be sure that is not used in calc
 * @param pX_m Pointer to the predicted output at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param pY_m Pointer to the propagetad sigma points array at (k|k-1) moment (i.e prediction in moment k based on states in (k-1))
 * @param sigmaIdx Sigma point index.
 */
void Hy2(tMatrix* pu, tMatrix* pX_m, tMatrix* pY_m, uint8_t sigmaIdx) {
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
