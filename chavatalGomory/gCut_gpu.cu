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

listNeigh *createGPUlist(const listNeigh* list_t){
    size_t size_list = sizeof(listNeigh) +
                        sizeof(TList)*(list_t->nList) +
                        sizeof(TPosList)*(list_t->nPos) ;

    listNeigh* h_list_gpu = (listNeigh*) malloc(size_list);
	memcpy(h_list_gpu,list_t, size_list);
	listNeigh* d_list;
	gpuMalloc((void**)&d_list, size_list);
	//printf("malloc ok\n");
	//getchar();
	gpuMemset(d_list,0,size_list);
	//printf("menset ok\n");
	//getchar();
	h_list_gpu->list_n = (TList*)(d_list + 1);
	h_list_gpu->pos = (TPosList*)(h_list_gpu->list_n + h_list_gpu->nList);

	gpuMemcpy(d_list, h_list_gpu, size_list, cudaMemcpyHostToDevice);
	free(h_list_gpu);
	return d_list;
}

