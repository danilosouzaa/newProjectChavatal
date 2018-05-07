#ifndef PREPARE_CPU_H
#define  PREPARE_CPU_H

//#include "cut_gpu.h"
//#include "configGpu.h"
#include "solutionGpu.h"
#include <sys/time.h>

EXTERN_C_BEGIN

void runCPUR1(Cut_gpu *h_cut, solutionGpu *h_solution, int precision);

Cut_gpu* initial_runCPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int maxDenominator, int precision);

void runCPUR2(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int *setConstraint, int precision, int maxDenominator, int nRuns, double timeGPU);

Cut_gpu* second_phase_runCPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int szR, double timeGPU);

EXTERN_C_END
#endif
