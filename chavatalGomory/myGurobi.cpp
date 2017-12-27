extern "C"{
    #include <stdio.h>
    #include <string.h>
    #include "cut_gpu.h"
    #include "Instance.h"
    #include "lp.h"

}



LinearProgramPtr geraLP(const char *fileName, Cut_gpu *cut_const, Instance *inst)
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
    for(i=0;i<number_variables;i++){
        name[i] = (char*)malloc(sizeof(char)*number_variables);
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
            c_variables[cont] = inst->Coef_obj[i];
            lb[cont]=0;
            ub[cont]=1;
            if(inst->type_variables[i]==1){
                integer[cont]=1;
            }else{
                integer[cont]=0;
            }
            strcpy(name[cont], inst->name_variables[i].name);
            cont++;
        //}
         indexes[i]=i;
    }

    lp_add_cols(lp,cont,c_variables,lb,ub,integer,name);


    for(i=0; i<inst->number_constraints; i++)
    {
        n_right[i] = inst->rhs[i];
    }

    char nome[20];
    for(i=0; i<inst->number_constraints; i++)
    {
        strcpy(nome, inst->name_constraints[i].name);
        cont=0;
        for(j=0;j<inst->number_variables;j++){
            if(inst->Coef_contraints[j + i*inst->number_variables] != 0){
                c_variables[cont] = inst->Coef_contraints[j + i*inst->number_variables];
                strcpy(name[cont], inst->name_variables[j].name);
                index_temp[cont] = indexes[j];
                cont++;

            }
        }
        if(inst->signal[i]==1){
            lp_add_row(lp,cont,index_temp,c_variables,nome,'E',n_right[i]);
        }else if(inst->signal[i]==4){
            lp_add_row(lp,cont,index_temp,c_variables,nome,'G',n_right[i]);

        }else if(inst->signal[i]==5){
            lp_add_row(lp,cont,index_temp,c_variables,nome,'L',n_right[i]);
        }else{
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
    return lp;

}



void runSolver(const char *fileName, int sizeFixSoft, float time){


}
