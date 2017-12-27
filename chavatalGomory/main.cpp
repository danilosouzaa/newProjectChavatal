extern "C"{
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
    int maxContraints = atoi(argv[3]);
    int nRuns = atoi(argv[4]);
    int maxDenominator = atoi(argv[5]);
    int type = atoi(argv[6]);
    printf("%d  \n",type);
    char fileName[50]	= "situation/";
    sprintf(fileName,"%s%s",fileName,nome);
    ccg = readFile(fileName,p);
    solutionGpu *sol = allocationStructSolution(ccg, maxContraints,nRuns);
    #ifdef __NVCC__
         printf("GPU\n");
         initial_runGPU(ccg, sol, maxContraints,nRuns,maxDenominator,p,type);
    #else
        printf("CPU\n");
        printf("Number Contraints: %d\n",ccg->numberConstrains);
    #endif // __NVCC__
    Instance* inst;
    inst = readLP("Danilo_teste.lp");
    printf("%d %d\n",inst->number_variables, inst->number_constraints);
    LinearProgramPtr lp = geraLP("Danilo_teste2.lp", ccg, inst);
    lp_free(&lp);
    printf("Hello world!\n");
    free(ccg);
    return 0;
}
