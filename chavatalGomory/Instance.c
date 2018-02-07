#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Instance.h"
#include "lp.h"

/*    int *ub_variables;
    int *lb_variables;
    short int *type_variables;*/

Instance *allocationStructInstance(int number_Variables, int number_Constraints)
{
    size_t size_inst = sizeof(Instance) +
                       sizeof(double)*number_Variables +
                       sizeof(short int)*number_Constraints +
                       sizeof(double)*number_Constraints*number_Variables+
                       sizeof(double)*number_Constraints +
                       sizeof(names)*number_Variables+
                       sizeof(names)*number_Constraints+
                       sizeof(double)*number_Variables+
                       sizeof(double)*number_Variables+
                       sizeof(short int)*number_Variables;

    Instance *inst = (Instance*)malloc(size_inst);
    assert(inst!=NULL);
    memset(inst,0,size_inst);
    inst->Coef_obj = (double*)(inst+1);
    inst->signal = (short int*)(inst->Coef_obj + number_Variables);
    inst->Coef_contraints = (double*)(inst->signal+ number_Constraints );
    inst->rhs = (double*)(inst->Coef_contraints + (number_Constraints*number_Variables));
    inst->name_variables = (names*)(inst->rhs + number_Constraints);
    inst->name_constraints = (names*)(inst->name_variables + number_Variables );
    inst->ub_variables = (double*)(inst->name_constraints + number_Constraints);
    inst->lb_variables = (double*)(inst->ub_variables+ number_Variables);
    inst->type_variables = (short int*)(inst->lb_variables+number_Variables);

    inst->number_constraints = number_Constraints;
    inst->number_variables= number_Variables;
    return inst;

}

int search_variables(Instance *inst, char *name)
{
    int i;
    for(i = 0; i < inst->number_variables; i++)
    {
        if(strcmp(inst->name_variables[i].name, name)==0)
        {
            return i;
        }
    }
    return -1;
}

Instance *readLP(char *fileName)
{
    FILE *arq;
    //arq = fopen(fileName,"r");
    char linha[255];
    int cont = 0, atv=0, cont2 = 0;

    double Coef,rhs, s=1;
    char *var;
    char signal[2];
    int nRes = 0, nCoef = 0;

    char *p,*p2;
    /*
    if(arq==NULL)
    {
        printf("Erro ao abrir o arquivo\n");
    }
    else
    {
        while(!feof(arq))
        {
            fscanf(arq,"%s",linha);
            //fscanf(arq,"%s\n",linha);
            //getchar();
            //getchar();
            if(strcmp(linha,"End")==0)
            {
                atv = 0;
            }
            if(strcmp(linha,"Generals")==0){
                cont--;
            }
            if(atv==1)
            {
                cont++;
            }
            p = strchr(linha,':');
            if(p!=NULL)
            {
                cont2++;
            }
            if(strcmp(linha,"Binaries")==0)
            {
                atv = 1;
            }
        }
        //printf("Numero de variaveis: %d\n", cont);
        //printf("Numero de restrições: %d\n", cont2);
        //getchar();
    }
    fclose(arq);
*/


    LinearProgramPtr lp;
    lp = lp_create();
    lp_read(lp,fileName);

    cont = lp_cols(lp);
    cont2 = lp_rows(lp);
    Instance *inst = allocationStructInstance(cont,cont2);
    printf("%d %d\n",inst->number_variables, inst->number_constraints);
    fflush(stdin);
//    getchar();
    int i;
    char* n_tes;
    for(i=0;i<cont;i++){

        n_tes = lp_col_name(lp,i,n_tes);
        strcpy(inst->name_variables[i].name, n_tes);
        nRes = lp_col_lb(lp,i);
        inst->lb_variables[i] = nRes;
        nRes = lp_col_ub(lp,i);
        inst->ub_variables[i] = nRes;
        if(lp_is_binary(lp,i)==1){
            inst->type_variables[i] = 1;
        }else if(lp_is_integer(lp,i)==1){
            inst->type_variables[i] = 2;
        }else{
            inst->type_variables[i] = 3;
        }
       // printf("%s %f %f %d\n", inst->name_variables[i].name, inst->lb_variables[i],inst->ub_variables[i],inst->type_variables[i] );
    }
    printf("%d %d\n",inst->number_variables, inst->number_constraints);
    for(i=0;i<cont2;i++){
        n_tes = lp_row_name(lp,i,n_tes);
        strcpy(inst->name_constraints[i].name, n_tes);

    }

    lp_free(&lp);
    int aux;
    arq = fopen(fileName,"r");
    if(arq==NULL)
    {
        printf("Erro ao abrir o arquivo\n");
    }
    else
    {
        nRes = 0;
        while(!feof(arq))
        {
            if(p==NULL)
            {
                fscanf(arq,"%s",linha);
                p2 = strchr(linha,':');
            }
            else
            {

                p=NULL;
            }

            aux = 0;

            if((strcmp(linha,"Minimize")==0)||(strcmp(linha,"Maxmize")==0)){
                while(aux==0){
                     fscanf(arq,"%s",linha);
                     if(strcmp(linha,"Subject")==0){
                        break;
                     }
                     if(strcmp(linha,"+")==0)
                    {
                        s = 1;
                        fscanf(arq,"%s",linha);
                        //printf("Positivo\n");
                    }
                    if(strcmp(linha,"-")==0)
                    {
                        s = -1;
                        fscanf(arq,"%s",linha);
                    }
                    if(( atof(linha) == 0 )&&(strcmp(linha,"0")!=0))
                    {
                        //printf("Não tem Coeficiente\n");
                        //getchar();
                        Coef = 1 * s;
                        var = linha;

                        i = search_variables(inst,var);
                        //printf("coef %f var %s position %d\n",Coef, var,i);
                        if(i!=-1){
                            inst->Coef_obj[i] = Coef;
                        }//getchar();
                    }else
                    {
                        Coef = atof(linha)*s;
                        fscanf(arq,"%s",linha);
                        var = linha;
                        i = search_variables(inst,var);
                        if(i != -1){
                            inst->Coef_obj[i] = Coef;
                        }
                        //getchar();
                    }


                }
            }


            if(p2!=NULL)
                fscanf(arq,"%s",linha);
            s = 1;
            while((p==NULL)&&(p2!=NULL)&&(strcmp(linha,"Bounds")!=0))
            {
                if((strcmp(linha,"=")!=0)&&(strcmp(linha,">")!=0)&&(strcmp(linha,"<")!=0)&&(strcmp(linha,">=")!=0)&&(strcmp(linha,"<=")!=0) )
                {
                    if(strcmp(linha,"+")==0)
                    {
                        s = 1;
                        fscanf(arq,"%s",linha);
                        //printf("Positivo\n");
                    }
                    if(strcmp(linha,"-")==0)
                    {
                        s = -1;
                        fscanf(arq,"%s",linha);
                    }
                    if( atof(linha) == 0 )
                    {
                        //printf("Não tem Coeficiente\n");
                        //getchar();
                        Coef = 1 * s;
                        var = linha;
                        i = search_variables(inst,var);
                        if(i != -1){

                        //printf("coef %f var %s position %d\n",Coef, var,i);
                        inst->Coef_contraints[i + nRes*inst->number_variables] = Coef;
                        }
                        //getchar();
                    }
                    else
                    {
                        Coef = atof(linha)*s;
                        fscanf(arq,"%s",linha);
                        var = linha;
                        i = search_variables(inst,var);
                       // printf("coef %f var %s position %d\n",Coef, var,i);
                        if(i != -1){


                            inst->Coef_contraints[i + nRes*inst->number_variables] = Coef;
                        }
                        //getchar();
                    }
                }
                else
                {


                    strcpy(signal,linha);
                    fscanf(arq,"%s",linha);
                    rhs = atof(linha);
                    if(nRes<inst->number_constraints)
                        inst->rhs[nRes] = rhs;
                    if(strcmp(signal,"=")==0){
                        inst->signal[nRes] = 1;
                    }

                    if(strcmp(signal,">")==0){
                        inst->signal[nRes] = 2;
                    }

                    if(strcmp(signal,"<")==0){
                        inst->signal[nRes] = 3;
                    }

                    if(strcmp(signal,">=")==0){
                        inst->signal[nRes] = 4;
                    }

                    if(strcmp(signal,"<=")==0){
                        inst->signal[nRes] = 5;
                    }
                    //printf("signal %s %d rhs %f\n",signal, inst->signal[nRes], rhs);
                    nRes++;
                    //getchar();
                }
                fscanf(arq,"%s",linha);
                p = strchr(linha,':');

            }


        }
//        printf("%s", inst->name_variables[0].name);
    }

    fclose(arq);
    return inst;
}
