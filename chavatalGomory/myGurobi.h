#ifndef MYGUROBI_H
#define  MYGUROBI_H

LinearProgramPtr geraLP(const char *fileName, Cut_gpu *cut_const, Instance *inst);

#endif
