#include <errno.h>
#include <libopencm3/cm3/vector.h>
#include <libopencm3/stm32/usart.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>

#include "clock.h"
#include "gpio.h"
#include "ukfLib/cfg/ukfCfg.h"
#include "usart.h"

#define UKF_TEST_EPS (1e-3)

extern tUkfMatrix UkfMatrixCfg;

void ukf_test(void);

int main(void) {
    clock_setup();
    gpio_setup();
    usart_setup();

    printf("App STARTED\n\n");

    //UKF test start here
    ukf_test();

    printf("\nApp DONE\n");

    while (1);

    return 0;
}

/******************************************************************************************************************************************************************************************************\
 ***  FUNCTION:
 ***      void ukf_test(void)
 *** 
 ***  DESCRIPTION:
 ***      Initialize and test UKF C implementation against expected result. Filter is tested in the loop from 15 steps. 
 ***      Total root square error is accumulated in the same loop for each state in order to show deviation from reference matlab solution.      
 ***            
 ***  PARAMETERS:
 ***      Type               Name              Range              Description
 ***      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ***      void
 ***  RETURNS:
 ***      void
 ***
 ***  SETTINGS:
 ***
\******************************************************************************************************************************************************************************************************/
void ukf_test(void) {
    uint8_t tfInitCfg = 0;
    tUKF ukfIo;
    uint32_t simLoop;

    //UKF filter measurement input(data log is generated in matlab and used for UKF simulation for 15 iteration)
    static const float yt[2][15] =
        {{0, 16.085992708563385, 12.714829185978214, 14.528500994457660, 19.105561355310275, 23.252820029388918, 29.282949862903255, 36.270058819651275, 44.244884173240955, 47.394243121124411, 55.988905459180458, 61.667450941562109, 68.624980301613647, 76.337963872393104, 82.611325690835159},  //y1 test
        {0,  16.750821420874981, 14.277640835870006, 16.320754051600520, 20.560460303503849, 24.827446289454556, 31.290961393448615, 36.853553457560210, 42.157283183453522, 49.382835230961490, 57.516319669684677, 65.664496283509095, 71.428712755732704, 79.241720894223079, 84.902760328915676}};

    //UKF filter expected system states calculated with matlab script for 15 iterations
    static const float x_exp[15][4] =
        {/*          x1                   x2                   x3                   x4*/
         {4.901482729572258,  4.576939885855807,  49.990342921246459, 49.958134463327802},
         {10.103304943868373, 9.409135720815829,  50.226544716205318, 49.750795004242228},
         {15.132069573131298, 14.138974122835807, 50.429540890147599, 49.191128327737864},
         {20.322823824348411, 19.096919763991380, 50.836860772439010, 49.189580207886742},
         {24.940146120267713, 23.399758647105461, 49.577386595072561, 47.383813382660449},
         {30.021901202161214, 27.882145089050120, 49.977568320123794, 46.551744562547626},
         {34.844137036519108, 32.753891693435087, 49.474006205027358, 47.190693993547214},
         {39.048783329419251, 38.499098203031146, 47.606199725375902, 50.001113730363919},
         {43.883085498256158, 42.383331307689538, 47.209657232695072, 46.747611757031784},
         {49.479941190207498, 47.255980559687778, 49.911944505272395, 47.887236233284476},
         {55.928745858553086, 51.180270357916882, 53.472542964944132, 46.510558543249353},
         {61.636426955126616, 55.275415649334157, 54.052126522632797, 45.262815392265203},
         {67.755622369652016, 59.602096868661732, 55.881393486598796, 45.326766104509289},
         {73.045763444967164, 63.838187852739992, 54.782159791340007, 44.291415099856643},
         {80.489525793047093, 66.908477563332085, 58.973616985147245, 42.638148924845950}};

    //UKF initialization: CFG
    tfInitCfg = ukf_init(&ukfIo, &UkfMatrixCfg);

    if (tfInitCfg == 0) {
        float err[4] = {0, 0, 0, 0};
        float absErrAccum[4] = {0, 0, 0, 0};

        //UKF simulation CFG0: BEGIN
        //printf("ukf.m | ukf.c | diff \n");
        for (simLoop = 1; simLoop < 15; simLoop++) {
            float* const py_cfg = ukfIo.input.y.val;

            //UKF:CFG0 apply/load system measurements in working array for current iteration.
            py_cfg[0] = yt[0][simLoop];
            py_cfg[1] = yt[1][simLoop];

            //UKF:CFG0 periodic task call
            (void)ukf_step(&ukfIo);

            err[0] = fabs(ukfIo.update.x.val[0] - x_exp[simLoop - 1][0]);
            err[1] = fabs(ukfIo.update.x.val[1] - x_exp[simLoop - 1][1]);
            err[2] = fabs(ukfIo.update.x.val[2] - x_exp[simLoop - 1][2]);
            err[3] = fabs(ukfIo.update.x.val[3] - x_exp[simLoop - 1][3]);

            /*
            printf("% -3.5f % -3.5f % -3.14f\n", x_exp[simLoop - 1][0], ukfIo.update.x.val[0], err[0]);
            printf("% -3.5f % -3.5f % -3.14f\n", x_exp[simLoop - 1][1], ukfIo.update.x.val[1], err[1]);
            printf("% -3.5f % -3.5f % -3.14f\n", x_exp[simLoop - 1][2], ukfIo.update.x.val[2], err[2]);
            printf("% -3.5f % -3.5f % -3.14f\n", x_exp[simLoop - 1][3], ukfIo.update.x.val[3], err[3]);
            */

            //accumulate the differennce between reference matlab implementation and results from C code execution
            absErrAccum[0] += err[0];
            absErrAccum[1] += err[1];
            absErrAccum[2] += err[2];
            absErrAccum[3] += err[3];
        }

        printf("Accumulated error between ukf.m and ukf.c \n");
        if (fabs(absErrAccum[0]) > UKF_TEST_EPS) { 
            printf("ERROR: Accumulated error absErrAccum[0] is too big: %.6e > %.6e\n", absErrAccum[0], UKF_TEST_EPS); 
        } else {
            printf("1. SUCCESS! %.6e < %.6e\n", absErrAccum[0], UKF_TEST_EPS); 
        } 
        if (fabs(absErrAccum[1]) > UKF_TEST_EPS) { 
            printf("ERROR: Accumulated error absErrAccum[1] is too big: %.6e > %.6e\n", absErrAccum[1], UKF_TEST_EPS); 
        } else {
            printf("2. SUCCESS! %.6e < %.6e\n", absErrAccum[1], UKF_TEST_EPS); 
        }
        if (fabs(absErrAccum[2]) > UKF_TEST_EPS) { 
            printf("ERROR: Accumulated error absErrAccum[2] is too big: %.6e > %.6e\n", absErrAccum[2], UKF_TEST_EPS); 
        } else {
            printf("3. SUCCESS! %.6e < %.6e\n", absErrAccum[2], UKF_TEST_EPS); 
        }
        if (fabs(absErrAccum[3]) > UKF_TEST_EPS) { 
            printf("ERROR: Accumulated error absErrAccum[3] is too big: %.6e > %.6e\n", absErrAccum[3], UKF_TEST_EPS); 
        } else {
            printf("4. SUCCESS! %.6e < %.6e\n", absErrAccum[3], UKF_TEST_EPS); 
        }

    } else {
        printf("initialization fail\n");
    }
}
