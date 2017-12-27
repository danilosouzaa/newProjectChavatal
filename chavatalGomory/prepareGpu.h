#ifndef PREPARE_GPU_H
#define  PREPARE_GPU_H

//#include "cut_gpu.h"
//#include "configGpu.h"
#include "solutionGpu.h"
#include "gpulib/types.h"




//void createCuts(Cut_gpu *h_cut, solutionGpu *h_solution,CutCG *ccg, int cont);
int verifyGpu();

void setGpuThread(int nGpu);

void initial_runGPU(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int nRuns, int maxDenominator, int precision, int type);




#endif
