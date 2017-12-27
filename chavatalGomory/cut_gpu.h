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

Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables);

Cut_gpu *readFile(char *fileName, int precision);

Cut_gpu* createGPUcut(const Cut_gpu* h_cut, int nVariables, int nConstrains);

EXTERN_C_END

#endif // CUT_GPU_H_
