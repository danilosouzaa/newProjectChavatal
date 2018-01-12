/*
 * gCut.cu
 *
 *  Created on: 29/03/2017
 *      Author: danilo
 */

#include "gCut_gpu.cuh"

Cut_gpu* createGPUcut(const Cut_gpu* h_cut, int nVariables, int nConstrains){
	//printf("begin createGPUInstance \n");
	int cont = h_cut->cont;
	size_t size_cut = sizeof(Cut_gpu) +
                        sizeof(TCoefficients)*(h_cut->cont) +
                        sizeof(TElements)*(h_cut->cont) +
                        sizeof(TElementsConstraints)*(nConstrains+1) +
                        sizeof(TRightSide)*(nConstrains) +
                        sizeof(TXAsterisc)*(nVariables)+
                        sizeof(TTypeConstraints)*(nConstrains);

	Cut_gpu* h_cut_gpu = (Cut_gpu*) malloc(size_cut);
	memcpy(h_cut_gpu,h_cut, size_cut);
	Cut_gpu* d_cut;
	gpuMalloc((void**)&d_cut, size_cut);
	//printf("malloc ok\n");
	//getchar();
	gpuMemset(d_cut,0,size_cut);
	//printf("menset ok\n");
	//getchar();
    h_cut_gpu->Coefficients = (TCoefficients*)(d_cut+1);
	h_cut_gpu->Elements = (TElements*)(h_cut_gpu->Coefficients + cont);
	h_cut_gpu->ElementsConstraints = (TElementsConstraints*)(h_cut_gpu->Elements + cont);
	h_cut_gpu->rightSide = (TRightSide*)(h_cut_gpu->ElementsConstraints + (nConstrains+1));
	h_cut_gpu->xAsterisc = (TXAsterisc*)(h_cut_gpu->rightSide + (nConstrains));
	h_cut_gpu->typeConstraints = (TTypeConstraints*)(h_cut_gpu->xAsterisc + (nVariables));
	h_cut_gpu->numberVariables = nVariables;
	h_cut_gpu->numberConstrains = nConstrains;
	h_cut_gpu->cont = cont;
	gpuMemcpy(d_cut, h_cut_gpu, size_cut, cudaMemcpyHostToDevice);
	free(h_cut_gpu);
	return d_cut;


}

