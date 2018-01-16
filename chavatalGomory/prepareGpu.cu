#include "gpulib/gpu.cuh"
//#include "gCut_gpu.cuh"
#include "gSolutionGpu.cuh"


extern "C" {
#include "prepareGpu.h"

}


void setGpuThread(int nGpu)
{
    gpuSetDevice(nGpu);
    int n;
    gpuGetDevice(&n);
    printf("gpu number %d\n", n);
}

int verifyGpu()
{
    int deviceCount = 0;
    //Commands for verify use correct of GPU
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
    if (error_id != cudaSuccess)
    {
        //printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
        //printf("Result = FAIL\n");
        return -1;
        //exit(1);
    }
    if(deviceCount == 0)
    {
        //printf("No GPU found :(");
        exit(1);
        return -1;
    }
    else
    {
        //printf("Found %d GPUs!\n", deviceCount);
        gpuSetDevice(0);
        //printf("GPU 0 initialized!\n");
        return deviceCount;
    }

}

Cut_gpu* initial_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int type)
{
    int deviceCuda;
    deviceCuda = verifyGpu();
    Cut_gpu* out_h_cut;
    if(deviceCuda > 0)
    {
        int i;//,j, nCons = h_cut->numberConstrains;
        size_t size_solution_r1 =  sizeof(solutionGpu) +
                                   sizeof(TSMult)*(nRuns) +
                                   sizeof(TSConst)*(nRuns) +
                                   sizeof(TSPAux)*(nRuns);

        size_t size_cut = sizeof(Cut_gpu) +
                          sizeof(TCoefficients)*(h_cut->cont) +
                          sizeof(TElements)*(h_cut->cont) +
                          sizeof(TElementsConstraints)*(h_cut->numberConstrains+1) +
                          sizeof(TRightSide)*(h_cut->numberConstrains) +
                          sizeof(TXAsterisc)*(h_cut->numberVariables) +
                          sizeof(TTypeConstraints)*(h_cut->numberConstrains);

        solutionGpu *h_solution_r1 = allocationStructSolution1(h_cut,nRuns);
        solutionGpu *d_solution_r1 = createGPUsolution1(h_solution_r1, h_cut,nRuns);
        Cut_gpu *d_cut = createGPUcut(h_cut, h_cut->numberVariables, h_cut->numberConstrains);
        curandState_t *states;
        cudaMalloc((void**)&states, (nRuns)*sizeof(curandState_t));
        unsigned int *h_seed = (unsigned int*)malloc(sizeof(unsigned int)*(nRuns));
        unsigned int *d_seed;
        srand(time(NULL));
        for(i=0; i<(nRuns); i++)
        {
            h_seed[i] = rand()%100000;
        }
        gpuMalloc((void*)&d_seed, sizeof(unsigned int)*(nRuns));
        gpuMemcpy(d_seed, h_seed, sizeof(unsigned int)*(nRuns), cudaMemcpyHostToDevice);

        int nT = nRuns/10;//nCons/10;
        int nB = 10;
        if(type==1)
            runGPUR1<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision);
        else
            runGPUR1_aleatory<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision, maxDenominator);

        gpuDeviceSynchronize();

        gpuMemcpy(h_solution_r1, d_solution_r1, size_solution_r1, cudaMemcpyDeviceToHost);
        h_solution_r1->SMult = (TSMult*)(h_solution_r1 + 1);
        h_solution_r1->SConst= (TSConst*)(h_solution_r1->SMult + (nRuns));
        h_solution_r1->SPAux = (TSPAux*)(h_solution_r1->SConst + (nRuns));


        gpuMemcpy(h_cut, d_cut, size_cut, cudaMemcpyDeviceToHost);
        h_cut->Coefficients = (TCoefficients*)(h_cut + 1);
        h_cut->Elements = (TElements*)(h_cut->Coefficients + (h_cut->cont));
        h_cut->ElementsConstraints = (TElementsConstraints*)(h_cut->Elements + (h_cut->cont));
        h_cut->rightSide = (TRightSide*)(h_cut->ElementsConstraints+ (h_cut->numberConstrains+1));
        h_cut->xAsterisc = (TXAsterisc*)(h_cut->rightSide + (h_cut->numberConstrains));
        h_cut->typeConstraints = (TTypeConstraints*)(h_cut->xAsterisc+ (h_cut->numberVariables));

        gpuFree(d_solution_r1);
        gpuFree(d_cut);
        gpuFree(d_seed);
        gpuFree(states);
        free(h_seed);

        int cont=0;

        //getchar();

        for(i=0; i<nRuns; i++)
        {
            if(h_solution_r1->SConst[i]!=-1)
            {
                //printf("%d %d /%d \n", h_solution_r1->SConst[i], h_solution_r1->SMult[i], h_solution_r1->SPAux[i]);
                cont++;
            }
        }

        if(cont>0)
        {
            printf("Number cuts generated in the phase 1: %d\n", cont);
            out_h_cut = createCutsOfPhaseOne(h_cut, cut_aux, h_solution_r1, cont,precision,nRuns);
            printf("DEPOIS!");
        }
        else
        {
            printf("No cuts generate\n");
            free(h_solution_r1);
            //free(h_cut);
            return h_cut;
        }

        free(h_solution_r1);
        free(h_cut);
    }

    return out_h_cut;

}

Cut_gpu* second_phase_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision)
{
    int deviceCuda;
    deviceCuda = verifyGpu();
    int *consR1;
    int *consNR1;
    int *nElemR1;

    Cut_gpu* out_cut_gpu;

    int n_r = 0, n_nr = 0, i;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==RES_R1)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            n_r++;
        }
        else
        {
            n_nr++;
        }
    }
    consR1 = (int*)malloc(sizeof(int)*n_r);
    nElemR1 = (int*)malloc(sizeof(int)*n_r);
    consNR1 = (int*)malloc(sizeof(int)*n_nr);

    n_r = 0;
    n_nr = 0;

    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==RES_R1)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            consR1[n_r] = i;
            nElemR1[n_r] = h_cut->ElementsConstraints[i+1] - h_cut->ElementsConstraints[i];
            n_r++;
        }
        else
        {
            consNR1[n_nr]=i;
            n_nr++;
        }
    }
    bubble_sort(nElemR1,consR1,n_r);
    int *Similar = returnOrdConstrainsNR(h_cut);
    float *folga = returnFolga(h_cut);

    solutionGpu *h_solution_r2 = allocationStructSolution2(h_cut,numberMaxConst,nRuns);
    int *setConstraint = (int*)malloc(sizeof(int)*numberMaxConst*nRuns);
    calcSetConstraint(setConstraint, numberMaxConst, h_cut->numberConstrains, consR1, consNR1, n_r, n_nr, Similar, folga,  nRuns );
    if(deviceCuda>0)
    {
        solutionGpu *d_solution;
        Cut_gpu *d_cut;
        int *d_setConstraint;

        int i, j, nB,nT;
        nB = 10;
        nT = nRuns/nB;

        size_t size_solution =  sizeof(solutionGpu) +
                                sizeof(TSMult)*(nRuns*4) +
                                sizeof(TSConst)*(numberMaxConst*nRuns) +
                                sizeof(TSPAux)*(nRuns);


        size_t size_cut = sizeof(Cut_gpu) +
                          sizeof(TCoefficients)*(h_cut->cont) +
                          sizeof(TElements)*(h_cut->cont) +
                          sizeof(TElementsConstraints)*(h_cut->numberConstrains+1) +
                          sizeof(TRightSide)*(h_cut->numberConstrains) +
                          sizeof(TXAsterisc)*(h_cut->numberVariables) +
                          sizeof(TTypeConstraints)*(h_cut->numberConstrains);

        d_solution = createGPUsolution2(h_solution_r2,h_cut,numberMaxConst,nRuns);
        d_cut = createGPUcut(h_cut,h_cut->numberVariables,h_cut->numberConstrains);

        curandState_t *states;
        cudaMalloc((void**)&states, (nT*nB)*sizeof(curandState_t));

        unsigned int *h_seed = (unsigned int*)malloc(sizeof(unsigned int)*(nT*nB));
        unsigned int *d_seed;
        srand(time(NULL));
        for(i=0; i<(nT*nB); i++)
        {
            h_seed[i] = rand()%100000;
        }
        gpuMalloc((void*)&d_seed, sizeof(unsigned int)*(nT*nB));
        gpuMemcpy(d_seed, h_seed, sizeof(unsigned int)*(nT*nB), cudaMemcpyHostToDevice);

        gpuMalloc((void*)&d_setConstraint, sizeof(int)*numberMaxConst*nRuns);
        gpuMemcpy(d_setConstraint, setConstraint, sizeof(int)*numberMaxConst*nRuns, cudaMemcpyHostToDevice);
        runGPUR2<<<nB,nT>>>(d_cut, d_solution, d_seed, states, numberMaxConst,d_setConstraint,nT,precision,maxDenominator);
        gpuDeviceSynchronize();
        gpuMemcpy(h_solution_r2, d_solution, size_solution, cudaMemcpyDeviceToHost);

        h_solution_r2->SMult = (TSMult*)(h_solution_r2+1);
        h_solution_r2->SConst= (TSConst*)(h_solution_r2->SMult + (nRuns*4));
        h_solution_r2->SPAux = (TSPAux*)(h_solution_r2->SConst + (numberMaxConst*nRuns));

        gpuMemcpy(h_cut, d_cut, size_cut, cudaMemcpyDeviceToHost);
        h_cut->Coefficients = (TCoefficients*)(h_cut+1);
        h_cut->Elements = (TElements*)(h_cut->Coefficients + h_cut->cont);
        h_cut->ElementsConstraints = (TElementsConstraints*)(h_cut->Elements + h_cut->cont);
        h_cut->rightSide = (TRightSide*)(h_cut->ElementsConstraints + (h_cut->numberConstrains+1));
        h_cut->xAsterisc = (TXAsterisc*)(h_cut->rightSide + (h_cut->numberConstrains));
        h_cut->typeConstraints = (TTypeConstraints*)(h_cut->xAsterisc + (h_cut->numberVariables));

        free(h_seed);
        gpuFree(states);
        gpuFree(d_setConstraint);
        gpuFree(d_cut);
        gpuFree(d_solution);
        gpuFree(d_seed);
        int cont=0;
        //printf("Number constraints: %d\n", h_cut->numberConstrains);
        for(i=0; i<nT; i++)
        {
            for(j=0; j<nB; j++)
            {
                if(h_solution_r2->SConst[0 + i*numberMaxConst + j*numberMaxConst*nT]!=-1)
                {
                    //printf("%d %d %d\n ",h_solution->SSize[i],h_solution->SPos[i],h_solution->SPAux[i]);
                    //printf("u1: %d / %d \t\t u2: %d / %d\n", h_solution->SMult[i], h_solution->SMult[i + 5*h_cut->numberConstrains], h_solution->SMult[i + 10*h_cut->numberConstrains], h_solution->SMult[i + 15*h_cut->numberConstrains]);
                    cont++;
                }
            }
        }
        if(cont>0)
        {
            printf("Number of Cuts in the second phase:%d\n",cont);
            out_cut_gpu = createCutsOfPhaseTwo(h_cut,cut_aux,h_solution_r2,numberMaxConst,cont,precision,nRuns,nT,nB);
        }
        else
        {
            free(consR1);
            free(consNR1);
            free(Similar);
            free(folga);
            free(nElemR1);
            free(setConstraint);
            free(h_solution_r2);

            return h_cut;

        }


    }
    free(consR1);
    free(consNR1);
    free(Similar);
    free(folga);
    free(nElemR1);
    free(setConstraint);
    free(h_solution_r2);
    free(h_cut);
    return out_cut_gpu;

}


