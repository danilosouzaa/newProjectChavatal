extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
//#include "cut_gpu.h"
#include "prepareGpu.h"
#include "Instance.h"
#include "lp.h"

}



#include "myGurobi.h"
//Primeiro Parametro: nome do arquivo
//Segundo Parametro: precisão de x
//Terceiro Parametro: número máximo de Restrição por blocos
//Quarto Parametro: número de execuções a serem realizadas
//Quinto Parametro: Número do maior denominador
//sexto parametro: Tipo do rank 1, aleatório ou nao.
//Setimo parametro: Size nr1
int main(int argc, char *argv[])
{
    if(argc<8)
    {
        printf("Argumentos faltando\n");
        return 0;
    }
    char *nome = argv[1];
    int p = atoi(argv[2]);
    Cut_gpu *ccg;
    Cut_gpu_aux *ccg_aux;
    int maxContraints = atoi(argv[3]);
    int nRuns = atoi(argv[4]);
    int maxDenominator = atoi(argv[5]);
    int type = atoi(argv[6]);
    int szR = atoi(argv[7]);

    if(szR>maxContraints){
        printf("Param erro: size R \n");
        return 0;
    }
    char fileName[50]	= "situation/";
    char nameInstance[50];
    int n1,n2;
    printf("%s\n",nome);
    sprintf(fileName,"%s%s",fileName,nome);
    int contr1 = 0,contr2 =0, n_cuts = 0 ;
    Instance* inst;

    struct timeval stop, start;
    gettimeofday(&start, NULL);
//do stuff


    FILE *arq;
    arq = fopen(fileName,"r");
    fscanf(arq,"%s \n %d %d %d\n",nameInstance,&n1,&n2,&contr1);
    fclose(arq);
    char nameLP[50];
    sprintf(nameLP,"lp/%s.lp",nameInstance);
    ccg_aux = AllocationStructCutAux(n1,n2);
    ccg = readFile(fileName,p,ccg_aux);


    fflush(stdin);
    inst = readLP(nameLP);
    LinearProgramPtr lp = geraLP("Danilo_teste2.lp", inst);


    contr1 = 0;
    int x=0;

#ifdef __NVCC__
        printf("GPU\n");
        int nBlocks;
        int pos_R1 = 0;
        int nThreads;
        int nRepeat = 1;
        int nRuns_temp,i;
        float aux;
        returnDimension(&nBlocks, &nThreads, nRuns);

        n_cuts= ccg->numberConstrains;
        //printf("antes: %d\n",ccg->numberConstrains);
        ccg = initial_runGPU(ccg, ccg_aux, maxContraints,maxDenominator,p,type,nThreads,nBlocks);
        //printf("depois fase 1: %d\n",ccg->numberConstrains);
        if(n_cuts!=ccg->numberConstrains)
            ccg_aux = reallocCut(ccg,ccg_aux, &contr1);
        n_cuts = ccg->numberConstrains;

        if(nRuns < nBlocks*nThreads){
            nThreads = (nRuns/nBlocks) ;
        }else{
            aux = (float)(nRuns - nThreads*nBlocks)/(float)(nThreads*nBlocks);
            nRepeat += ceil(aux);
        }
        for(i = 1; i <= nRepeat;i++){
            if((nRepeat>1)&&(i == nRepeat)){
                nThreads = (nRuns - (nThreads*nBlocks)*(i-1))/nBlocks;
            }
            nRuns_temp = nThreads*nBlocks;
            ccg = second_phase_runGPU(ccg, ccg_aux, maxContraints,nRuns_temp,maxDenominator,p, nBlocks, nThreads, &pos_R1, szR);
                //printf("depois fase 2: %d - %d\n",ccg->numberConstrains, pos_R1);
            //if(n_cuts!=ccg->numberConstrains)
            //    ccg_aux = reallocCutR2(ccg,ccg_aux,&contr2);
        }
        ccg_aux = reallocCutR2(ccg,ccg_aux,&contr2);
        //ccg = phase_zeroHalf(ccg, ccg_aux,2);
#else
        printf("CPU\n");
        printf("Number Contraints: %d\n",ccg->numberConstrains);
#endif // __NVCC__

        gettimeofday(&stop, NULL);
        double secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
        printf("took %f\n", secs);

        lp = InsertCutsInLP(lp,ccg,ccg_aux,inst);
        lp_set_max_seconds(lp,10);
        int yy = lp_optimize_as_continuous(lp);
        //updateXAstherisc(ccg,ccg_aux,lp,p);

        //yy = lp_optimize(lp);
        //x++;

    lp_free(&lp);

    free(ccg_aux);
    free(ccg);

    //free(col);
    return 0;
}
