#ifndef MYGUROBI_H
#define  MYGUROBI_H

LinearProgramPtr geraLP(const char *fileName, Instance *inst);

int get_number_constraints(const char *fileName);

int get_number_variables(const char *fileName);

int searchConstLP(char *name,TNames *namesConst, int szNamesConst);

int returnIndex(char *name, TNames *namesVar, int szNamesVar);

LinearProgramPtr InsertCutsInLP(LinearProgramPtr lp,Cut_gpu *ccg, Cut_gpu_aux *ccg_aux, Instance *inst);

void updateXAstherisc(Cut_gpu *ccg,Cut_gpu_aux *ccg_aux, LinearProgramPtr lp,int precision);
#endif
