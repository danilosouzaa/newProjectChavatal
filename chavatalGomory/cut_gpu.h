#ifndef CUT_GPU_H_
#define CUT_GPU_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "gpulib/types.h"

EXTERN_C_BEGIN

typedef int TCoefficients;
typedef int TElements;
typedef int TElementsConstraints;
typedef int TRightSide;
typedef int TXAsterisc;
typedef int TTypeConstraints;
typedef int TInterval;

typedef int TList;
typedef int TPosList;


typedef struct{
    char name[255];
}TNames;

typedef struct {
    int numberVariables;
    int numberConstrains;
    int cont;
    TCoefficients *Coefficients;
    TElements *Elements;
    TElementsConstraints *ElementsConstraints;
    TRightSide *rightSide;
    TXAsterisc *xAsterisc;
    TTypeConstraints *typeConstraints;

}Cut_gpu;


typedef struct{
    int numberVariables;
    int numberConstrains;
    TInterval *intervalMin;
    TInterval *intervalMax;
    TNames *nameElements;
    TNames *nameConstraints;
}Cut_gpu_aux;


typedef struct {
    int nList;
    int nPos;
    TList *list_n;
    TPosList *pos;
}listNeigh;



Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables);

Cut_gpu_aux *AllocationStructCutAux(int nConstrains, int nVariables);

listNeigh *AllocationListNeigh(int nConstrains, int nList);

Cut_gpu* readFile(char *fileName, int precision ,Cut_gpu_aux *cut_aux);

//Cut_gpu *readFile(char *fileName, int precision, Cut_gpu_aux* cut_aux);

Cut_gpu* createGPUcut(const Cut_gpu* h_cut, int nVariables, int nConstrains);

listNeigh *createGPUlist(const listNeigh* list_t, int nConstrains, int nList);

int returnIndVector(TNames *v,char *nome, int sz);

EXTERN_C_END

#endif // CUT_GPU_H_
