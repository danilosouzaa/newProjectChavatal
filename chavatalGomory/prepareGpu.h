#ifndef PREPARE_GPU_H
#define  PREPARE_GPU_H

//#include "cut_gpu.h"
//#include "configGpu.h"
#include "solutionGpu.h"
#include "gpulib/types.h"
#include <sys/time.h>


//void createCuts(Cut_gpu *h_cut, solutionGpu *h_solution,CutCG *ccg, int cont);
int verifyGpu();

void setGpuThread(int nGpu);

//Cut_gpu* initial_runGPU(Cut_gpu *h_cut,Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int type);

Cut_gpu* initial_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int maxDenominator, int precision, int type, int nThreads, int nBlocks);

//Cut_gpu* second_phase_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision);

Cut_gpu* second_phase_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int nB,int nT, int *pos_R1, int szR);




Cut_gpu* phase_zeroHalf(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux,int nConstraintsPerSet);

listNeigh *returnMatrixNeighborhood (Cut_gpu *h_cut);

int contPar(Cut_gpu* h_cut);

void *shuffle_Set(int *vec, int nSetConstrains, int n);


void returnDimension(int *nB, int *nT, int nRuns);

void fillParImpar(int *vPar,int *vImpar, Cut_gpu *h_cut);
#endif
