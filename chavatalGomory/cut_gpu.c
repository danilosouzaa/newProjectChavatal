#include "cut_gpu.h"

Cut_gpu *AllocationStructCut(int cont, int nConstrains, int nVariables)
{
    size_t size_cut = sizeof(Cut_gpu) +
                      sizeof(TCoefficients)*(cont) +
                      sizeof(TElements)*(cont) +
                      sizeof(TElementsConstraints)*(nConstrains+1) +
                      sizeof(TRightSide)*(nConstrains) +
                      sizeof(TXAsterisc)*(nVariables)+
                      sizeof(TTypeConstraints)*(nConstrains);
    Cut_gpu *cut = (Cut_gpu*)malloc(size_cut);
    assert(cut!=NULL);
    memset(cut,0,size_cut);
    cut->Coefficients = (TCoefficients*)(cut+1);
    cut->Elements = (TElements*)(cut->Coefficients + cont);
    cut->ElementsConstraints = (TElementsConstraints*)(cut->Elements + cont);
    cut->rightSide = (TRightSide*)(cut->ElementsConstraints + (nConstrains+1));
    cut->xAsterisc = (TXAsterisc*)(cut->rightSide + (nConstrains));
    cut->typeConstraints = (TTypeConstraints*)(cut->xAsterisc + (nVariables));
    cut->numberVariables = nVariables;
    cut->numberConstrains = nConstrains;
    cut->cont = cont;
    return cut;
}



Cut_gpu *readFile(char *fileName, int precision)
{
    FILE *arq;
    char nameInstance[50];
    Cut_gpu *ccg;
    char sense ='+';
    arq = fopen(fileName,"r");
    int n1,n2,i,cont;
    float coef;
    int rhs;
    if(arq==NULL)
    {
        printf("Erro ao abrir o arquivo\n");
    }
    else
    {
        int pos = 0, elem = 0, temp = 0;
        fscanf(arq,"%s \n %d %d %d\n",nameInstance,&n1,&n2,&cont);
        printf("ts: %s, t1: %d, t2: %d cont: %d\n", nameInstance,n1,n2,cont);
        getchar();
        ccg = AllocationStructCut(cont,n1,n2);
        ccg->ElementsConstraints[0] = 0;
//	ccg->numberConstrains = n1;
//	ccg->numberVariables = n2;
        for(i = 0; i < n1; i++)
        {
            while(sense=='+')
            {
                fscanf(arq,"%f x%d %c",&coef, &elem, &sense);
                ccg->Coefficients[pos] = coef;
                ccg->Elements[pos] = elem;
                pos++;
                //printf(" %f %d %c", coef,cont,sense);
            }

            fscanf(arq," = %d",&rhs);
            ccg->rightSide[i] = rhs;
            ccg->ElementsConstraints[i+1] = pos;
            sense = '+';
        }
        for(i=0; i<n2; i++)
        {
            fscanf(arq,"\nx*%d = %f\n",&elem,&coef);
            ccg->xAsterisc[elem] = precision*coef;
        }
        for(i=0; i<n1; i++)
        {
            fscanf(arq,"\nr%d %d\n",&elem,&pos);
            ccg->typeConstraints[i] = pos;
            printf("%d\n",pos);
        }
    }

    fclose(arq);
    return ccg;
}
