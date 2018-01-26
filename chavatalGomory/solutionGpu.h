/*
 * solutionGpu.h
 *
 *  Created on: 31/03/2017
 *      Author: danilo
*/

#ifndef SOLUTION_GPU_H_
#define SOLUTION_GPU_H_


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "cut_gpu.h"
#include "gpulib/types.h"
#include "lp.h"

EXTERN_C_BEGIN

//typedef int TSCoefficients;
//typedef int TSRightSide;
//typedef int TSViolation;
typedef int TSMult;
typedef int TSConst;
typedef int TSPAux;


typedef struct {
    TSMult *SMult;
    TSConst *SConst;
    TSPAux *SPAux;
} solutionGpu;

solutionGpu* allocationStructSolution2(Cut_gpu *c, int numberMaxConst, int nRuns);

solutionGpu* allocationStructSolution1(Cut_gpu *c, int nRuns);

solutionGpu* createGPUsolution2(solutionGpu* h_solution, Cut_gpu* h_cut,int numberMaxConst, int nRuns);

solutionGpu* createGPUsolution1(solutionGpu* h_solution, Cut_gpu* h_cut, int nRuns);




Cut_gpu* createCutsOfPhaseOne(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, solutionGpu *h_solution, int nCuts, int precision, int nRuns);

Cut_gpu* createCutsOfPhaseTwo(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, solutionGpu *h_solution, int numberMaxConst, int nCuts, int precision, int nRuns, int nThreads,int nBlocks);

int CutP_maxDivisorCommonVector(int coefs[], int nElem);

int CutP_maxDivisorCommonRec(int m, int n);

Cut_gpu_aux* reallocCut(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux, int *cont);

Cut_gpu_aux* reallocCutR2(Cut_gpu *h_cut,Cut_gpu_aux *h_cut_aux, int *cont);

void bubble_sort(int *vetor,int *pos, int n);

int* returnOrdConstrainsNR(Cut_gpu *cut);

float* returnFolga(Cut_gpu *cut);

//void calcSetConstraint (int *setConstraint, int numberMaxConst,int numberConstrains, int *resR1, int *resNR1, int sizeR1, int sizeNR1, int *Similar, float *Folga, int nRuns);
void calcSetConstraint (int *setConstraint, int *pos_R1, int numberMaxConst,int numberConstrains, int *resR1, int *resNR1, int sizeR1, int sizeNR1, int *Similar, float *Folga, int nRuns );


EXTERN_C_END

#endif

