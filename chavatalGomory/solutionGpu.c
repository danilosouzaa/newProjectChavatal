/*
 * solutionGpu.c
 *
 *  Created on: 29/03/2017
 *      Author: danilo
*/

#include "solutionGpu.h"




solutionGpu* allocationStructSolution(Cut_gpu *c, int numberMaxConst, int nRuns)
{
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns*4) +
                            sizeof(TSConst)*(numberMaxConst*nRuns) +
                            sizeof(TSPAux)*(nRuns);
    solutionGpu *sol;
    sol = (solutionGpu*)malloc(size_solution);
    assert(sol!=NULL);
    memset(sol,0,size_solution);
    sol->SMult = (TSMult*)(sol+1);
    sol->SConst= (TSConst*)(sol->SMult + (nRuns*4));
    sol->SPAux = (TSPAux*)(sol->SConst + (numberMaxConst*nRuns));
    return sol;
}

solutionGpu* allocationStructSolution1(Cut_gpu *c, int nRuns)  //for phase 1
{
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns) +
                            sizeof(TSConst)*(nRuns) +
                            sizeof(TSPAux)*(nRuns);
    solutionGpu *sol;
    sol = (solutionGpu*)malloc(size_solution);
    assert(sol!=NULL);
    memset(sol,0,size_solution);
    sol->SMult = (TSMult*)(sol+1);
    sol->SConst= (TSConst*)(sol->SMult + (nRuns));
    sol->SPAux = (TSPAux*)(sol->SConst + (nRuns));
    return sol;
}

