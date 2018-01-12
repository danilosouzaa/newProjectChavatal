#ifndef INSTANCE_H
#define  INSTANCE_H



typedef struct{
    char name[255];
}names;

//Signal:
// 1 =
// 2 >
// 3 <
// 4 >=
// 5 <=

typedef struct{
    int number_variables; //ok
    int number_constraints;//ok
    double *Coef_obj; //ok
    short int *signal;//ok
    double *Coef_contraints;//ok
    double *rhs;//ok
    names *name_variables;//ok
    names *name_constraints;//ok
    double *ub_variables;
    double *lb_variables;
    short int *type_variables; //(1 for binaries,and 2 for integer, 3 real)
}Instance;


Instance *allocationStructInstance(int number_Variables, int number_Constraints);

int search_variables(Instance *inst, char *name);

Instance *readLP(char *fileName);


#endif
