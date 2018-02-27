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


solutionGpu* createGPUsolution2(solutionGpu* h_solution, Cut_gpu* h_cut,int numberMaxConst, int nRuns)
{

    size_t size_solution =  sizeof(solutionGpu) +
                            sizeof(TSMult)*(nRuns*4) +
                            sizeof(TSConst)*(numberMaxConst*nRuns) +
                            sizeof(TSPAux)*(nRuns);

    solutionGpu *d_sol;
    gpuMalloc((void**)&d_sol, size_solution);
    gpuMemset(d_sol,0,size_solution);
    h_solution->SMult = (TSMult*)(d_sol+1);
    h_solution->SConst= (TSConst*)(h_solution->SMult + (nRuns*4));
    h_solution->SPAux = (TSPAux*)(h_solution->SConst + (numberMaxConst*nRuns));
    gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
    return d_sol;
}


__global__ void runGPUR1(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int nThreads, int precision)
{


    int term = threadIdx.x + blockIdx.x*nThreads;
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
    int term = threadIdx.x + blockIdx.x*nThreads;
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



__global__ void runGPUR2(Cut_gpu *d_cut, solutionGpu *d_solution, unsigned int *seed, curandState_t* states, int numberMaxConst, int setConstraint[],int nThreads, int precision, int maxDenominator, int nRuns)
{
    int term = threadIdx.x + blockIdx.x*nThreads;
    if(term<nRuns)
    {
       // printf("%d: %d %d %d %d\n",term, setConstraint[term*numberMaxConst + 0],setConstraint[term*numberMaxConst + 1],setConstraint[term*numberMaxConst + 2],setConstraint[term*numberMaxConst + 3]);
        int mult_1, mult_2, rest_a,rest_b, i, j, el, rhs1, rhs2, value_tes, violation = 0, aux, n1_best = -1, n2_best = -1, d1_best = -1, qnt_1, d2_best=-1;//, cont=0;
        curand_init(seed[term],term,0,&states[term]);
        int Numerator[20];
        int Denominator[20];
        int *Coef = (int*)malloc(sizeof(int)*(d_cut->numberVariables));
        int *Coef2 = (int*)malloc(sizeof(int)*(d_cut->numberVariables));
        for(i=0; i<20; i++)
        {
            Denominator[i]= curand(&states[term])%maxDenominator + 2;
            Numerator[i] = curand(&states[term])%(Denominator[i]-1);
        }
        for(mult_1=0; mult_1<20; mult_1++)
        {
            memset(Coef,0,sizeof(int)*d_cut->numberVariables);
            rhs1 = 0;
            for(rest_a = 0; rest_a< numberMaxConst; rest_a++)
            {
                for(i=d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_a] ]; i<d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_a] + 1]; i++)
                {

                    el = d_cut->Elements[i];
                    Coef[el] += d_cut->Coefficients[i] * Numerator[mult_1];
                }
                rhs1 += d_cut->rightSide[ setConstraint[term*numberMaxConst+rest_a] ] * Numerator[mult_1];
                for(mult_2 = 0; mult_2<20; mult_2++)
                {
                    memset(Coef2,0,sizeof(int)*d_cut->numberVariables);
                    value_tes = 0;
                    rhs2 = 0;
                    for(rest_b = rest_a + 1; rest_b<numberMaxConst; rest_b++)
                    {
                        for(j=d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_b] ]; j<d_cut->ElementsConstraints[ setConstraint[term*numberMaxConst + rest_b] + 1]; j++)
                        {
                            el = d_cut->Elements[j];
                            Coef2[el] += d_cut->Coefficients[j] * Numerator[mult_2];
                        }
                        rhs2 += d_cut->rightSide[ setConstraint[term*numberMaxConst + rest_b] ]* Numerator[mult_2];
                    }
                    for(j=0; j<d_cut->numberVariables; j++)
                    {
                        aux = Coef[j]<0 ? Coef[j]/Denominator[mult_1] - 1 : Coef[j]/Denominator[mult_1];
                        value_tes += aux*d_cut->xAsterisc[j];
                        aux = Coef2[j]<0 ? Coef2[j]/Denominator[mult_2] - 1 : Coef2[j]/Denominator[mult_2];
                        value_tes += aux*d_cut->xAsterisc[j];
                    }
                    aux = rhs1<0 ? rhs1/Denominator[mult_1]-1 : rhs1/Denominator[mult_1];
                    aux +=  rhs2<0 ? rhs2/Denominator[mult_2]-1 : rhs2/Denominator[mult_2];


                    if((value_tes>aux*precision)&&(value_tes-(aux*precision)>violation))
                    {
                        violation = value_tes-(aux*precision);
                        if(violation>precision)
                        {
                            printf("AQUIIII!!");
                            for(i=0; i<numberMaxConst; i++)
                            {
                                printf("%d ",setConstraint[term*numberMaxConst + i]);//CPU ja vai ter
                            }
                            printf("\n");
                        }
                        n1_best = Numerator[mult_1];
                        d1_best = Denominator[mult_1];
                        n2_best = Numerator[mult_2];
                        d2_best = Denominator[mult_2];
                        qnt_1 = rest_a;
                    }


                }
            }

        }
        __syncthreads();

        if(violation>0)
        {
            for(i=0; i<numberMaxConst; i++)
            {
                d_solution->SConst[i + threadIdx.x*numberMaxConst + blockIdx.x*numberMaxConst*nThreads] = setConstraint[term*numberMaxConst + i];//CPU ja vai ter
            }

            d_solution->SPAux[threadIdx.x + blockIdx.x*nThreads] = qnt_1;
            d_solution->SMult[0 + threadIdx.x*4 + blockIdx.x*4*nThreads] = n1_best;
            d_solution->SMult[1 + threadIdx.x*4 + blockIdx.x*4*nThreads] = d1_best;
            d_solution->SMult[2 + threadIdx.x*4 + blockIdx.x*4*nThreads] = n2_best;
            d_solution->SMult[3 + threadIdx.x*4 + blockIdx.x*4*nThreads] = d2_best;

        }
        else
        {
            for(i=0; i<numberMaxConst; i++)
            {
                d_solution->SConst[i + threadIdx.x*numberMaxConst + blockIdx.x*numberMaxConst*nThreads] = -1;
            }
            d_solution->SPAux[threadIdx.x + blockIdx.x*nThreads] = 0;
            d_solution->SMult[0 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
            d_solution->SMult[1 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
            d_solution->SMult[2 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
            d_solution->SMult[3 + threadIdx.x*4 + blockIdx.x*4*nThreads] = -1;
        }



        free(Coef);
        free(Coef2);
        __syncthreads();
    }

}
