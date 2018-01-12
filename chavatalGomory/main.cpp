extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

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
int main(int argc, char *argv[])
{
    if(argc<7)
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
    char fileName[50]	= "situation/";
    char nameInstance[50];
    int n1,n2;

    sprintf(fileName,"%s%s",fileName,nome);
    int contr1 = 0,contr2 =0, n_cuts = 0 ;
    Instance* inst;
    fflush(stdin);
    inst = readLP("Danilo_teste.lp");
    LinearProgramPtr lp = geraLP("Danilo_teste2.lp", inst);

    FILE *arq;
    arq = fopen(fileName,"r");
    fscanf(arq,"%s \n %d %d %d\n",nameInstance,&n1,&n2,&contr1);
    fclose(arq);
    ccg_aux = AllocationStructCutAux(n1,n2);
    ccg = readFile(fileName,p,ccg_aux);
    contr1 = 0;
    int x=0;
    while(x<5)
    {
        //printf("%s\n", ccg_aux->nameConstraints[1].name);
        //solutionGpu *sol = allocationStructSolution(ccg, maxContraints,nRuns);
#ifdef __NVCC__
        printf("GPU\n");
        n_cuts= ccg->numberConstrains;
        printf("antes: %d\n",ccg->numberConstrains);
        ccg = initial_runGPU(ccg, ccg_aux, maxContraints,nRuns,maxDenominator,p,type);
        printf("depois fase 1: %d\n",ccg->numberConstrains);
        if(n_cuts!=ccg->numberConstrains)
            ccg_aux = reallocCut(ccg,ccg_aux, &contr1);
        n_cuts = ccg->numberConstrains;
        ccg = second_phase_runGPU(ccg, ccg_aux, maxContraints,nRuns,maxDenominator,p);
        printf("depois fase 2: %d\n",ccg->numberConstrains);
        if(n_cuts!=ccg->numberConstrains)
            ccg_aux = reallocCutR2(ccg,ccg_aux,&contr2);
#else
        printf("CPU\n");
        printf("Number Contraints: %d\n",ccg->numberConstrains);
#endif // __NVCC__
        lp = InsertCutsInLP(lp,ccg,ccg_aux,inst);

        int yy = lp_optimize_as_continuous(lp);
        updateXAstherisc(ccg,ccg_aux,lp,p);
        lp_set_max_seconds(lp,10);
        yy = lp_optimize(lp);
        x++;
    }
    lp_free(&lp);

    //free(ccg_aux);
    free(ccg);

    //free(col);
    return 0;
}
