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

void initial_runGPU(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int nRuns, int maxDenominator, int precision, int type)
{
    int deviceCuda;
    deviceCuda = verifyGpu();
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
        Cut_gpu *d_cut = createGPUcut(h_cut,h_cut->numberVariables, h_cut->numberConstrains);
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

        int nT = 50;//nCons/10;
        int nB = 20;
        if(type==1)
            runGPUR1<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision);
        else
            runGPUR1_aleatory<<<nB,nT>>>(d_cut, d_solution_r1, d_seed, states, nT, precision, maxDenominator);
        gpuDeviceSynchronize();

        gpuMemcpy(h_solution_r1, d_solution_r1, size_solution_r1, cudaMemcpyDeviceToHost);
        h_solution_r1->SMult = (TSMult*)(h_solution_r1+1);
        h_solution_r1->SConst= (TSConst*)(h_solution_r1->SMult + (nRuns));
        h_solution_r1->SPAux = (TSPAux*)(h_solution_r1->SConst + (nRuns));
        gpuMemcpy(h_cut, d_cut, size_cut, cudaMemcpyDeviceToHost);

        h_cut->Coefficients = (TCoefficients*)(h_cut+1);
        h_cut->Elements = (TElements*)(h_cut->Coefficients + h_cut->cont);
        h_cut->ElementsConstraints = (TElementsConstraints*)(h_cut->Elements + h_cut->cont);
        h_cut->rightSide = (TRightSide*)(h_cut->ElementsConstraints + (h_cut->numberConstrains+1));
        h_cut->xAsterisc = (TXAsterisc*)(h_cut->rightSide + (h_cut->numberConstrains));
        h_cut->typeConstraints = (TTypeConstraints*)(h_cut->xAsterisc + (h_cut->numberVariables));


        gpuFree(d_solution_r1);
        gpuFree(d_cut);
        gpuFree(d_seed);
        gpuFree(states);
        free(h_seed);
        int cont=0;
        //printf("Number constraints: %d\n number runs: %d\n", h_cut->numberConstrains,nRuns);
        //getchar();
        nRuns =nB*nT;
        for(i=0; i<nRuns; i++)
        {
                if(h_solution_r1->SConst[i]!=-1)
                {
                    //printf("%d %d /%d \n", h_solution_r1->SConst[i], h_solution_r1->SMult[i], h_solution_r1->SPAux[i]);
                    cont++;
                }
        }
        printf("Number cuts: %d\n", cont);
        /*CutCG *ccg_r1;
        if(cont>0)
        {
            printf("Preview Size of constraints: %d\n", VDbl_size(ccg->rowrhs));
            ccg_r1 = createCutsR1(h_cut, h_solution_r1,ccg,cont);
            printf("new Size of constraints: %d\n", VDbl_size(ccg_r1->rowrhs));
            printf("\nPhase 1: \t Num cuts gpu - %d\n", cont);
            int a1;
            VecInt *aux,*aux2, *timeInst;
            aux = VInt_create();
            aux2 = VInt_create();
            timeInst = VInt_create();


            printf("Initial of phase 2\n");
            for(i = 0 ; i<VDbl_size(ccg_r1->rowrhs); i++)
            {
                a1 = VInt_get(ccg_r1->rowtype,i);

                if((a1==RES_RR)||(a1==RES_R1))
                {
                    VInt_pushBack(aux,i);
                    VInt_pushBack(timeInst,VInt_get(ccg_r1->cutnelem,i));
                }
                else
                {
                    VInt_pushBack(aux2,i);
                }
            }

            int *resR1;
            int *timeR1;
            int *resNR1;
            resR1 = VInt_getPtr(aux);
            resNR1 = VInt_getPtr(aux2);
            timeR1 = VInt_getPtr(timeInst);
            bubble_sort(timeR1,resR1, VInt_size(aux));

            Cut_gpu *cutG;
            solutionGpu *solG;
            int cont_2 = 0;


            for( int r = 0 ; r < VDbl_size(ccg_r1->rowrhs) ; r++)
                cont_2 += VDbl_size(ccg_r1->rowCoef[r]);

            cutG =  AllocationStructCut(cont_2,VDbl_size(CutCG_getrowrhs(ccg_r1)),VDbl_size(CutCG_getxfElemPP(ccg_r1)));
            fillStructCut_gpu(cutG,ccg_r1);
            int *Similar =  returnOrdConstrainsNR(cutG,ccg_r1, resNR1);
            float *folga = returnFolga(cutG,ccg_r1, resNR1);
            // showStructCut_GPU(cutG)

            solG = allocationStructSolution(cutG,numberMaxConst);

            //runGPU_New(cutG,h_solution,ccg, numberMaxConst, 0);
            runGPU_phase2(cutG,solG, ccg, numberMaxConst, resR1, resNR1, VInt_size(aux), VInt_size(aux2),Similar,folga);

            // CutCG_printCut(ccg);
            //getchar();

            free(cutG);
            free(resR1);
            free(timeR1);
            free(resNR1);
            free(Similar);
            free(folga);

            CutCG_free(&ccg_r1);
        }
        else
        {
            printf("nothing cuts generated\n");
        }
        free(h_solution_r1);
        */
    }



}
