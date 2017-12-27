/*
 * gSolution.cu
 *
 *  Created on: 31/03/2017
 *      Author: danilo
 */
#include "gSolutionGpu.cuh"




solutionGpu* createGPUsolution1(solutionGpu* h_solution, Cut_gpu* h_cut, int nRuns)
{

    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns) +
                            sizeof(TSConst)*(nRuns) +
                            sizeof(TSPAux)*(nRuns);

    solutionGpu *d_sol;
    gpuMalloc((void**)&d_sol, size_solution);
    gpuMemset(d_sol,0,size_solution);
    h_solution->SMult = (TSMult*)(d_sol+1);
    h_solution->SConst= (TSConst*)(h_solution->SMult + (nRuns));
    h_solution->SPAux = (TSPAux*)(h_solution->SConst + (nRuns));
    gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
    return d_sol;
}






solutionGpu* createGPUsolution(solutionGpu* h_solution, Cut_gpu* h_cut,int numberMaxConst)
{
    int nThreads = 1;
    int nBlocks=2;
    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nThreads*nBlocks*4) +
                            sizeof(TSConst)*(numberMaxConst*nThreads*nBlocks) +
                            sizeof(TSPAux)*(nThreads*nBlocks);

    solutionGpu *d_sol;
    gpuMalloc((void**)&d_sol, size_solution);
    gpuMemset(d_sol,0,size_solution);
    h_solution->SMult = (TSMult*)(d_sol+1);
    h_solution->SConst= (TSConst*)(h_solution->SMult + (nThreads*nBlocks*4));
    h_solution->SPAux = (TSPAux*)(h_solution->SConst + (numberMaxConst*nBlocks*nThreads));
    gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
    return d_sol;
}


__global__ void runGPUR1(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision)
{


#define term  threadIdx.x + blockIdx.x*nThreads
    __shared__ int *constraints;
    __shared__ int pos;
    curand_init(seed[term],term,0,&states[term]);

    int violation = 0,i,j;
    if(threadIdx.x == 0)
    {
        pos = 0;
        constraints = (int*)malloc(sizeof(int)*d_cut->numberConstrains);
        for(i=0; i<d_cut->numberConstrains; i++)
        {
            if(d_cut->typeConstraints[i] == RES_RR)
            {
                constraints[pos] = i;
                pos++;
            }
        }

    }
    __syncthreads();

    int res = constraints[threadIdx.x%pos];
    int *Coef = (int*)malloc(sizeof(int)*(d_cut->numberVariables));

    int n1=-1, d1=-1,el, rhs, aux,value_tes;
    int nBest=-1, dBest=-1, violation_best=0;
    for(j = d_cut->ElementsConstraints[ res ] ; j < d_cut->ElementsConstraints[ res +1 ]; j++)
    {
        d1 = d_cut->Coefficients[j];
        n1 = 1;
        while(n1<d1)
        {
            rhs = 0;
            violation = 0;
            value_tes = 0;
            for(i = d_cut->ElementsConstraints[ res ]; i<d_cut->ElementsConstraints[ res + 1 ]; i++)
            {
                el = d_cut->Elements[i];
                aux = d_cut->Coefficients[i] * n1;
                if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                {
                    aux = (aux/d1) -1;
                }
                else
                {
                    aux = aux/d1;
                }
                //aux = aux< 0 ? (aux/d1) - 1 : aux/d1;
                value_tes += aux*d_cut->xAsterisc[el];
            }
            rhs = d_cut->rightSide[ res ]* n1;
            if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
            {
                rhs = (rhs/d1) -1;
            }
            else
            {
                rhs = rhs/d1;
            }

            if(value_tes>rhs*precision)
            {
                violation = value_tes - (rhs*precision);
                if(violation>violation_best)
                {
                    violation_best = violation;
                    nBest=n1;
                    dBest=d1;
                }
            }
            n1++;
        }
    }

    if(violation_best!=0)
    {
        d_solution->SConst[term] = res;
        d_solution->SMult[term] = nBest;
        d_solution->SPAux[term] = dBest;
    }
    else
    {
        d_solution->SConst[term] = -1;
        d_solution->SMult[term] = -1;
        d_solution->SPAux[term] = -1;
    }

    free(Coef);
    if(threadIdx.x == 0)
    {
        free(constraints);
    }
}


__global__ void runGPUR1_aleatory(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision,int maxDenominator)
{

#define term  threadIdx.x + blockIdx.x*nThreads
    __shared__ int *constraints;
    __shared__ int pos;
    curand_init(seed[term],term,0,&states[term]);

    int violation = 0, cont = 0,i;
    if(threadIdx.x == 0)
    {
        pos = 0;
        constraints = (int*)malloc(sizeof(int)*d_cut->numberConstrains);
        for(i=0; i<d_cut->numberConstrains; i++)
        {
            if(d_cut->typeConstraints[i] == RES_RR)
            {
                constraints[pos] = i;
                pos++;
            }
        }

    }
    __syncthreads();

    int res = constraints[threadIdx.x%pos];
    int *Coef = (int*)malloc(sizeof(int)*(d_cut->numberVariables));
    cont = 0;
    int n1=-1, d1=-1,el, rhs, aux,value_tes;
    int nBest=-1, dBest=-1, violation_best=0;
    while((cont<20)&&(violation_best==0))
    {
        cont++;
        d1 = curand(&states[term])%maxDenominator + 2;
        n1 = 1;
        while(n1<d1)
        {
            rhs = 0;
            violation = 0;
            value_tes = 0;
            //printf("%d/%d\n",n1,d1);
            for(i = d_cut->ElementsConstraints[ res ]; i<d_cut->ElementsConstraints[ res + 1 ]; i++)
            {
                el = d_cut->Elements[i];
                aux = d_cut->Coefficients[i] * n1;
                if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                {
                    aux = (aux/d1) -1;
                }
                else
                {
                    aux = aux/d1;
                }
                value_tes += aux*d_cut->xAsterisc[el];
            }
            rhs = d_cut->rightSide[ res ]* n1;
            if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
            {
                rhs = (rhs/d1) -1;
            }
            else
            {
                rhs = rhs/d1;
            }

            if(value_tes>rhs*precision)
            {
                violation = value_tes - (rhs*precision);
                if(violation>violation_best)
                {
                    violation_best = violation;
                    nBest=n1;
                    dBest=d1;
                }
            }
            n1++;
        }
    }

    if(violation_best!=0)
    {
        d_solution->SConst[term] = res;
        d_solution->SMult[term] = nBest;
        d_solution->SPAux[term] = dBest;
    }
    else
    {
        d_solution->SConst[term] = -1;
        d_solution->SMult[term] = -1;
        d_solution->SPAux[term] = -1;
    }

    free(Coef);
    if(threadIdx.x == 0)
    {
        free(constraints);
    }
}
