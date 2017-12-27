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

solutionGpu* allocationStructSolution(Cut_gpu *c, int numberMaxConst, int nRuns);

solutionGpu* allocationStructSolution1(Cut_gpu *c, int nRuns);

solutionGpu* createGPUsolution(solutionGpu* h_solution, Cut_gpu* h_cut, int numberMaxConst);

solutionGpu* createGPUsolution1(solutionGpu* h_solution, Cut_gpu* h_cut, int nRuns);




EXTERN_C_END

#endif

