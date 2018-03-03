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

void shuffle_Set(int *vec, int nSetConstrains, int n)
{
    timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    int i, j, aux ;
    int *num_temp = (int*)malloc(sizeof(int)*nSetConstrains);
    int *vec_aux = (int*)malloc(sizeof(int)*nSetConstrains);
    aux  =  n/nSetConstrains;
    for(i = 0; i < aux ; i++)
    {

        for(j = 0 ; j<nSetConstrains; j++)
        {

            num_temp[j] = rand()%RAND_MAX;
            vec_aux[j] = vec[i*nSetConstrains + j];
        }
        bubble_sort(num_temp,vec_aux,nSetConstrains);
        for(j = 0 ; j<nSetConstrains; j++)
        {
            vec[i*nSetConstrains + j] = vec_aux[j];
        }
    }
    free(num_temp);
    free(vec_aux);
}

Cut_gpu* initial_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int maxDenominator, int precision, int type, int nThreads, int nBlocks)
{
    int deviceCuda;
    cudaSetDevice(0);
    deviceCuda = verifyGpu();
    Cut_gpu* out_h_cut;
    int nRuns;
    if(deviceCuda > 0)
    {
        int i, numberC = 0 ;//,j, nCons = h_cut->numberConstrains;

        for(i = 0; i<h_cut->numberConstrains; i++)
        {
            if(h_cut->typeConstraints[i] == RES_RR)
            {
                numberC++;
            }
        }
        //float auxD = ((float)numberC)/((float)nBlocks);
        int nT = numberC;//nCons/10;
        //int nT = ceil(auxD);//nCons/10;
        int nB = 1;
        //int nB = nBlocks;
        nRuns = nT*nB;
//        nRuns = 1000;
//        nB = 10;
//        nT = 100;
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

        if(type==1)
        {
            runGPUR1<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision);
        }
        else
        {
            runGPUR1_aleatory<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision, maxDenominator);
        }
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


void returnDimension(int *nB, int *nT, int nRuns)
{

    int blockSize;      // The launch configurator returned block size
    int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch
    int gridSize;
    int N = nRuns;

    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,runGPUR2, 0, N);
    *nB = minGridSize;
    *nT = blockSize;
}


Cut_gpu* second_phase_runGPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int nB,int nT, int *pos_R1, int szR)
{
    int deviceCuda;
    deviceCuda = verifyGpu();
    int *consR1;
    int *consNR1;
    int *nElemR1;
    cudaSetDevice(0);
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
            if(h_cut->typeConstraints[i]!=LPC_CGGPUR2)
            {
                consNR1[n_nr]=i;
                n_nr++;
            }
        }
    }
    bubble_sort(nElemR1,consR1,n_r);
    int *Similar = returnOrdConstrainsNR(h_cut);
    float *folga = returnFolga(h_cut);
//    printf("%d %d %d\n",nRuns, nT, nB);
//
//    getchar();
    solutionGpu *h_solution_r2 = allocationStructSolution2(h_cut,numberMaxConst,nRuns);
    int *setConstraint = (int*)malloc(sizeof(int)*numberMaxConst*nRuns);
    calcSetConstraint(setConstraint, pos_R1,numberMaxConst, h_cut->numberConstrains, consR1, consNR1, n_r, n_nr, Similar, folga,  nRuns, szR);
    /*int j;
     for(i=0;i<nRuns;i++){
         for(j=0;j<numberMaxConst;j++){
             printf("%d \t", setConstraint[i*numberMaxConst + j]);

         }
         printf("\n");

     }
    */

    shuffle_Set(setConstraint, numberMaxConst, numberMaxConst*nRuns);

    /*  printf("Depois do Shuffle\n");
       for(i=0;i<nRuns;i++){
          printf("%d:", i);
          for(j=0;j<numberMaxConst;j++){
              printf("%d \t", setConstraint[i*numberMaxConst + j]);

          }
          printf("\n");
      }
    */

    if(deviceCuda>0)
    {
        solutionGpu *d_solution;
        Cut_gpu *d_cut;
        int *d_setConstraint;

        int i, j;
//        if(blockSize*minGridSize < nRuns){
//            nRp = nRuns - blockSize*minGridSize;
//        }

        //nB = 10;
        //nT = nRuns/nB;

        // nB = minGridSize;
        // nT = blockSize;

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
        gpuMalloc((void**)&d_seed, sizeof(unsigned int)*(nT*nB));
        gpuMemcpy(d_seed, h_seed, sizeof(unsigned int)*(nT*nB), cudaMemcpyHostToDevice);

        gpuMalloc((void*)&d_setConstraint, sizeof(int)*(numberMaxConst*nRuns));
        gpuMemcpy(d_setConstraint, setConstraint, sizeof(int)*(numberMaxConst*nRuns), cudaMemcpyHostToDevice);

        runGPUR2<<<nB,nT>>>(d_cut, d_solution, d_seed, states, numberMaxConst, d_setConstraint, nT,precision,maxDenominator,nRuns);
	//cudaSetDevice(0);
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
            if(out_cut_gpu==NULL)
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

int contPar(Cut_gpu* h_cut)
{
    int cont = 0,i;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->rightSide[i]%2==0)
        {
            cont++;
        }
    }
    return cont;
}

Cut_gpu* phase_zeroHalf(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux,int nConstraintsPerSet)
{
    //char *matrixNeighborhood;
    int i,j;
    int szPar = contPar(h_cut);
    int szImpar = h_cut->numberConstrains -szPar;
    int *vPar = (int*)malloc(sizeof(int)*szPar);
    int *vImpar = (int*)malloc(sizeof(int)*(szImpar));
    //matrixNeighborhood = returnMatrixNeighborhood(h_cut);

    fillParImpar(vPar,vImpar,h_cut);
    int nBlocks, nThreads;
    nBlocks = 10;
    nThreads =  szPar/nBlocks;
    int deviceCuda;
    deviceCuda = verifyGpu();
    if(deviceCuda>0)
    {
        int *h_solution_r2 = (int*)malloc(sizeof(int)*nConstraintsPerSet*nBlocks*nThreads);

        free(h_solution_r2);
    }

    free(vImpar);
    free(vPar);
    //free(matrixNeighborhood);
    return h_cut;

}

void fillParImpar(int *vPar,int *vImpar, Cut_gpu *h_cut)
{
    int i, cP=0, cI = 0;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->rightSide[i]%2==0)
        {
            vPar[cP] = i;
            cP++;
        }
        else
        {
            vImpar[cI] = i;
            cI++;
        }
    }
}


listNeigh *returnMatrixNeighborhood (Cut_gpu *h_cut)
{
    char *matrixNeighborhood = (char*)malloc(sizeof(char)*h_cut->numberConstrains*h_cut->numberConstrains);
    int *m1 = (int*)malloc(sizeof(int)*h_cut->numberConstrains*h_cut->numberVariables);
    int i,j, k, el, cont_temp = 0;
    memset(m1,0, sizeof(int)*h_cut->numberConstrains*h_cut->numberVariables);
    memset(matrixNeighborhood,0, sizeof(char)*h_cut->numberConstrains*h_cut->numberConstrains);
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        for(j = h_cut->ElementsConstraints[i]; j<h_cut->ElementsConstraints[i+1]; j++)
        {
            el = h_cut->Elements[j];
            m1[el + i*h_cut->numberVariables] = h_cut->Coefficients[j];
        }
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        for(j=0; j<h_cut->numberConstrains; j++)
        {
            for(k = 0; k<h_cut->numberVariables; k++)
            {
                if((i!=j)&&( ((m1[k + i*h_cut->numberVariables]>0)&&(m1[k + j*h_cut->numberVariables]>0)) || ((m1[k + i*h_cut->numberVariables]<0)&&(m1[k + j*h_cut->numberVariables]<0)) ) )
                {
                    matrixNeighborhood[i+j*h_cut->numberConstrains] = 1;
                    cont_temp++;
                    break;
                }
            }
        }
    }
    listNeigh *list_t = AllocationListNeigh(h_cut->numberConstrains,cont_temp);
    //int *novaLista = (int*)malloc(sizeof(int)*cont_temp);
    //int *pos = (int*)malloc(sizeof(int)*h_cut->numberConstrains+1);
    cont_temp = 0;
    list_t->pos[0] = 1;
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        for(j=0; j<h_cut->numberConstrains; j++)
        {
            if(matrixNeighborhood[i+j*h_cut->numberConstrains] == 1)
            {
                list_t->list_n[cont_temp] = j;
                cont_temp++;
            }
        }
        list_t->pos[i+1] = cont_temp;
    }

    free(m1);

    free(matrixNeighborhood);
    return list_t;
}


