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


Cut_gpu_aux *AllocationStructCutAux(int nConstrains, int nVariables){
     size_t size_cut_aux = sizeof(Cut_gpu_aux) +
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TInterval)*(nConstrains)+
                          sizeof(TNames)*(nVariables)+
                          sizeof(TNames)*(nConstrains);


    Cut_gpu_aux *cut_aux = (Cut_gpu_aux*)malloc(size_cut_aux);
    assert(cut_aux!=NULL);
    memset(cut_aux,0,size_cut_aux);
    cut_aux->intervalMin = (TInterval*)(cut_aux + 1);
    cut_aux->intervalMax = (TInterval*)(cut_aux->intervalMin + (nConstrains));
    cut_aux->nameElements = (TNames*)(cut_aux->intervalMax + (nConstrains));
    cut_aux->nameConstraints = (TNames*)(cut_aux->nameElements+(nVariables));
    cut_aux->numberVariables = nVariables;
    cut_aux->numberConstrains = nConstrains;
    return cut_aux;

}

listNeigh *AllocationListNeigh(int nConstrains, int nList){
    size_t size_list = sizeof(listNeigh) +
                          sizeof(TList)*(nList)+
                          sizeof(TPosList)*(nConstrains+1);
    listNeigh *list_t = (listNeigh*)malloc(size_list);
    assert(list_t!=NULL);
    memset(list_t, 0 ,size_list);
    list_t->list_n = (TList*)(list_t + 1);
    list_t->pos = (TPosList*)(list_t->list_n + nList);
    list_t->nList = nList;
    list_t->nPos = nConstrains + 1;
    return list_t;
}

int returnIndVector(TNames *v,char *nome, int sz){
    int i;
    for(i=0;i<sz;i++){
        if(strcmp(v[i].name,nome)==0)
            return i;
    }
    return -1;
}



Cut_gpu* readFile(char *fileName, int precision,Cut_gpu_aux *cut_aux)
{
    FILE *arq;
    char nameInstance[50];
    Cut_gpu *ccg;
    //Cut_gpu_aux *cut_aux;
    char linha[255];
    char sense ='+';
    arq = fopen(fileName,"r");
    int n1,n2,i,cont, a1, a2;
    float coef;
    int rhs;
    if(arq==NULL)
    {
        printf("Erro ao abrir o arquivo\n");
    }
    else
    {
        int pos = 0, elem = 0, temp = 0, sig = 1;
        fscanf(arq,"%s \n %d %d %d\n",nameInstance,&n1,&n2,&cont);
        //printf("ts: %s, t1: %d, t2: %d cont: %d\n", nameInstance,n1,n2,cont);
        //getchar();
        TNames *name_temp;
        name_temp = (TNames*)malloc(sizeof(TNames)*cont);
        ccg = AllocationStructCut(cont,n1,n2);
        //cut_aux = AllocationStructCutAux(n1,n2);


        ccg->ElementsConstraints[0] = 0;
//	ccg->numberConstrains = n1;
//	ccg->numberVariables = n2;
        for(i = 0; i < n1; i++)
        {
            fscanf(arq,"%s: ",linha);

            strncpy(cut_aux->nameConstraints[i].name,linha,strlen(linha)-1);

            //getchar();
            fscanf(arq,"%s",linha);
            if((atof(linha)==0.0)&&(strcmp(linha,"0")!=0))
            {
                sig = -1;
                fscanf(arq,"%f",&coef);
                coef = sig*coef;
            }
            else
            {
                coef = atof(linha);
            }
            ccg->Coefficients[pos] = coef;


            fscanf(arq,"%s %c",linha,&sense);

            //linha será a variavel
            strcpy(name_temp[pos].name,linha);
            pos++;
            while((sense=='+')||(sense=='-'))
            {
                if(sense=='-')
                {
                    sig = -1;
                }
                else
                {
                    sig = 1 ;
                }
                fscanf(arq,"%f %s %c",&coef, linha, &sense);
                ccg->Coefficients[pos] = sig*coef;
                strcpy(name_temp[pos].name,linha);
                //ccg->Elements[pos] = elem;
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
            fscanf(arq,"\n %s = %f \n",linha, &coef);
            strcpy(cut_aux->nameElements[i].name,linha);
            ccg->xAsterisc[i] = precision*coef;
        }
        for(i=0; i<n1; i++)
        {
            fscanf(arq,"\n%s %d %d %d\n",linha,&pos, &a1, &a2);
            ccg->typeConstraints[i] = pos;
            cut_aux->intervalMin[i] = a1;
            cut_aux->intervalMax[i] = a2;
            //printf("%d %d %d \n",ccg->typeConstraints[i], ccg->intervalMax[i], ccg->intervalMin[i]);
        }
        for(i=0;i<cont;i++){
            temp = returnIndVector(cut_aux->nameElements,name_temp[i].name,n2);
            ccg->Elements[i] = temp;
          //  printf("temp: %d\n",temp);
        }
       //getchar();
        free(name_temp);
    }
    fclose(arq);

    //Cut_Master *ccg_master;


    return ccg;
}
