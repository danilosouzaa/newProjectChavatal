extern "C" {
#include <stdio.h>
#include <string.h>
#include "cut_gpu.h"
#include "Instance.h"
#include "lp.h"

}

int get_number_variables(const char *fileName)
{
    LinearProgramPtr lp;
    lp = lp_create();
    lp_read(lp,fileName);
    int n_var = lp_cols(lp);
    lp_free(&lp);
    return n_var;
}

int get_number_constraints(const char *fileName)
{
    LinearProgramPtr lp;
    lp = lp_create();
    lp_read(lp,fileName);
    int n_cons = lp_rows(lp);
    lp_free(&lp);
    return n_cons;


}



LinearProgramPtr geraLP(const char *fileName, Instance *inst)
{
    LinearProgramPtr lp;//
    lp = lp_create();
    int number_variables = inst->number_variables; //cut_const->numberVariables;
    char n[15];
    strcpy(n, fileName);
    double *c_variables;
    double *lb;
    double *ub;
    char *integer; //descobrir para que serve.
    char **name;
    int *indexes;
    double *n_right;
    int i,j;
    int *index_temp;

    integer = (char*)malloc(sizeof(char)*number_variables);

    c_variables = (double*)malloc(sizeof(double)*number_variables);
    lb = (double*)malloc(sizeof(double)*number_variables);
    ub = (double*)malloc(sizeof(double)*number_variables);
    name = (char**)malloc(sizeof(char*)*number_variables);
    for(i=0; i<number_variables; i++)
    {
        name[i] = (char*)malloc(sizeof(char)*255);
    }
    indexes = (int*)malloc(sizeof(int)*number_variables);
    index_temp = (int*)malloc(sizeof(int)*number_variables);
    //cof = (double*)malloc(sizeof(double)*inst->nJobs);
    // cof_temp = (double*)malloc(sizeof(double)*inst->mAgents);
    n_right = (double*)malloc(sizeof(double)*(inst->number_constraints));
    int cont = 0;
    for(i=0; i<number_variables; i++)
    {
        //if(inst->Coef_obj[i] != 0){
        //printf("%f\n",inst->Coef_obj[i]);
        c_variables[i] = inst->Coef_obj[i];
        lb[i]=inst->lb_variables[i];
        ub[i]=inst->ub_variables[i];
        if((inst->type_variables[i]==1)||(inst->type_variables[i]==2))
        {
            integer[i]=1;
        }
        else
        {
            integer[i]=0;
        }
        //Relaxação
        /*if(inst->type_variables[i]==1){
            integer[i]=0;
        }*/

        strcpy(name[i], inst->name_variables[i].name);
        //cont++;
        //printf("branco: %s\t", inst->name_variables[i].name);
        //getchar();
        //}
        indexes[i]=i;
    }
    //getchar();
    lp_add_cols(lp,number_variables,c_variables,lb,ub,integer,name);
    printf("OK");
    //getchar();
    for(i=0; i<inst->number_constraints; i++)
    {
        n_right[i] = inst->rhs[i];
    }

    char nome[20];
    for(i=0; i<inst->number_constraints; i++)
    {
        strcpy(nome, inst->name_constraints[i].name);
        cont=0;
        for(j=0; j<inst->number_variables; j++)
        {
            if(inst->Coef_contraints[j + i*inst->number_variables] != 0)
            {
                c_variables[cont] = inst->Coef_contraints[j + i*inst->number_variables];
                strcpy(name[cont], inst->name_variables[j].name);
                index_temp[cont] = indexes[j];
                cont++;

            }
        }
        if(inst->signal[i]==1)
        {
            lp_add_row(lp,cont,index_temp,c_variables,nome,'E',n_right[i]);
        }
        else if(inst->signal[i]==4)
        {
            lp_add_row(lp,cont,index_temp,c_variables,nome,'G',n_right[i]);

        }
        else if(inst->signal[i]==5)
        {
            lp_add_row(lp,cont,index_temp,c_variables,nome,'L',n_right[i]);
        }
        else
        {
            printf("No equalite, problem\n");
        }
    }





//    for(i=0; i<number_variables; i++)
//    {
//        name[i]=new char[15];
//    }
//
//
    for(i=0; i<number_variables; i++)
    {
        free(name[i]);
    }
    free(name);
    free(lb);
    free(ub);

    free(integer);
    free(c_variables);
    free(n_right);
    free(indexes);
    free(index_temp);

    lp_write_lp(lp,fileName);
    lp_set_max_seconds(lp,300);

    return lp;

}

int searchConstLP(char *name,TNames *namesConst, int szNamesConst)
{
    int i,aux;
    for(i=0; i<szNamesConst; i++)
    {
        aux = strcmp(namesConst[i].name,name);
        if(aux==0)
        {
            return 1;
        }
    }
    return 0;
}


int returnIndex(char *name, TNames *namesVar, int szNamesVar){
    int i,aux;
    for(i=0; i<szNamesVar; i++)
    {
        aux = strcmp(namesVar[i].name,name);
        if(aux==0)
        {
            return i;
        }
    }
    return -1;

}

LinearProgramPtr InsertCutsInLP(LinearProgramPtr lp,Cut_gpu *ccg, Cut_gpu_aux *ccg_aux, Instance *inst)
{
    TNames *namesVar = (TNames*)malloc(sizeof(TNames)*inst->number_variables);

    int i;
    int szConstLp = lp_rows(lp);
    TNames *namesConst = (TNames*)malloc(sizeof(TNames)*szConstLp);

    char nameConst[255];
    int tes;
    for(i=0; i<szConstLp; i++)
    {
        lp_row_name( lp, i, nameConst );
        strcpy(namesConst[i].name,nameConst);
    }
    for(i=0;i<inst->number_variables;i++){
        lp_col_name(lp,i,nameConst);
        strcpy(namesVar[i].name,nameConst);
    }
    int cont;
    double *c_variables = (double*)malloc(sizeof(double)*ccg->numberVariables);
    int *index_temp = (int*)malloc(sizeof(int)*ccg->numberVariables);
    int j, c_aux;
    int el, idx;
    for(i=0; i<ccg_aux->numberConstrains; i++)
    {   strcpy(nameConst,ccg_aux->nameConstraints[i].name);
        tes = searchConstLP(nameConst,namesConst,szConstLp);
        if(tes==0)
        {
            cont = ccg->ElementsConstraints[i+1] - ccg->ElementsConstraints[i];
            c_aux = 0;
            for(j = ccg->ElementsConstraints[i];j<ccg->ElementsConstraints[i+1];j++){
                el= ccg->Elements[j];
                c_variables[c_aux] = ccg->Coefficients[j];
                idx = returnIndex(ccg_aux->nameElements[el].name, namesVar, inst->number_variables);
                if(idx!=-1){
                    index_temp[c_aux] = idx;
                    c_aux++;
                }else{
                    printf("Element Not found");
                }
            }
            lp_add_row(lp,cont,index_temp,c_variables,nameConst,'L',ccg->rightSide[i]);
        }
    }

    free(c_variables);
    free(namesVar);
    free(namesConst);
    lp_write_lp(lp,"danilo_ccs.lp");

    return lp;
}

void updateXAstherisc(Cut_gpu *ccg,Cut_gpu_aux *ccg_aux, LinearProgramPtr lp, int precision)
{
    double *x_astherisc;

    int szVariables,i, idx;
    szVariables = lp_cols(lp);
    x_astherisc = lp_x(lp);
    char name[255];
    TNames *namesVariables = (TNames*)malloc(sizeof(TNames)*szVariables);
    for(i=0; i<szVariables; i++)
    {
        lp_col_name(lp,i,name);
        strcpy(namesVariables[i].name,name);
    }
    for(i=0;i<ccg->numberConstrains;i++){
        idx = returnIndex(ccg_aux->nameElements[i].name,namesVariables,szVariables);
        ccg->xAsterisc[i] = x_astherisc[idx]*precision;
    }

}


