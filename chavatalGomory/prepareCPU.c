#include "prepareCPU.h"

void runCPUR1(Cut_gpu *h_cut, solutionGpu *h_solution, int precision)
{
    int *constraints;
    int pos;
    int violation = 0, i, j;
    pos = 0;
    constraints = (int*)malloc(sizeof(int)*h_cut->numberConstrains);
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i] == RES_RR) //&& (d_cut->typeConstraints[i] == LPC_ENERGETIC))

        {
            constraints[pos] = i;
            pos++;

        }
    }
    int k;
    for(k = 0; k<pos; k++)
    {
        int res = constraints[k];
        int *Coef = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
        int n1=-1, d1=-1,el, rhs, aux,value_tes;
        int nBest=-1, dBest=-1, violation_best=0;
        for(j = h_cut->ElementsConstraints[ res ] ; j < h_cut->ElementsConstraints[ res +1 ]; j++)
        {
            d1 = h_cut->Coefficients[j];
            n1 = 1;
            while(n1<d1)
            {
                rhs = 0;
                violation = 0;
                value_tes = 0;
                for(i = h_cut->ElementsConstraints[ res ]; i<h_cut->ElementsConstraints[ res + 1 ]; i++)
                {
                    el = h_cut->Elements[i];
                    aux = h_cut->Coefficients[i] * n1;
                    if( ((aux>0&&d1<0)||(aux<0&&d1>0))&&(aux%d1!=0))
                    {
                        aux = (aux/d1) -1;
                    }
                    else
                    {
                        aux = aux/d1;
                    }
                    //aux = aux< 0 ? (aux/d1) - 1 : aux/d1;
                    value_tes += aux*h_cut->xAsterisc[el];
                }
                rhs = h_cut->rightSide[ res ]* n1;
                if( ((rhs>0&&d1<0)||(rhs<0&&d1>0))&&(rhs%d1!=0))
                {
                    rhs = (rhs/d1) -1;
                }
                else
                {
                    rhs = rhs/d1;
                }

                if(value_tes>rhs*precision)
                {
                    violation = value_tes - (rhs*precision);
                    if(violation>violation_best)
                    {
                        violation_best = violation;
                        nBest=n1;
                        dBest=d1;
                    }
                }
                n1++;
            }
        }

        if(violation_best!=0)
        {
            h_solution->SConst[k] = res;
            h_solution->SMult[k] = nBest;
            h_solution->SPAux[k] = dBest;
        }
        else
        {
            h_solution->SConst[k] = -1;
            h_solution->SMult[k] = -1;
            h_solution->SPAux[k] = -1;
        }
        free(Coef);
    }
    free(constraints);
}


Cut_gpu* initial_runCPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int maxDenominator, int precision)
{
    Cut_gpu* out_h_cut;
    int nRuns;
    int i, numberC = 0 ;
    for(i = 0; i<h_cut->numberConstrains; i++)
    {
        if(h_cut->typeConstraints[i] == RES_RR) //&& (d_cut->typeConstraints[i] == LPC_ENERGETIC))

        {
            numberC++;
        }
    }
    nRuns = numberC;
    solutionGpu *h_solution_r1 = allocationStructSolution1(h_cut,nRuns);
    runCPUR1(h_cut, h_solution_r1,  precision);
    int cont=0;
    for(i=0; i<nRuns; i++)
    {
        if(h_solution_r1->SConst[i]!=-1)
        {
            cont++;
        }
    }

    if(cont>0)
    {
        printf("Number cuts generated in the phase 1: %d\n", cont);
        //out_h_cut = createCutsStrongPhaseOne(h_cut, h_solution_r1, cont,precision,nRuns);
        out_h_cut = createCutsOfPhaseOne(h_cut,cut_aux, h_solution_r1, cont,precision,nRuns);

        free(h_solution_r1);
        free(h_cut);
        return out_h_cut;
    }
    else
    {
        printf("No cuts generate\n");
        free(h_solution_r1);
        return h_cut;
    }



}


void runCPUR2(Cut_gpu *h_cut, solutionGpu *h_solution, int numberMaxConst, int *setConstraint, int precision, int maxDenominator, int nRuns,double timeGPU)
{
    int mult_1, mult_2, rest_a,rest_b, i, j, el, rhs1, rhs2, value_tes, violation = 0, aux, n1_best = -1, n2_best = -1, d1_best = -1, qnt_1, d2_best=-1;//, cont=0;
    double timeCurrent = omp_get_wtime();
    int Numerator[20];
    int Denominator[20];
    int k;
    int *Coef = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
    int *Coef2 = (int*)malloc(sizeof(int)*(h_cut->numberVariables));
//    printf("TIME: %lf\n",timeGPU);
//    timeGPU = ((double) timeGPU - (omp_get_wtime() - timeCurrent ));
//    if(timeGPU<1)
//    {
//        for(k=0; k<nRuns; k++)
//        {
//            for(i=0; i<numberMaxConst; i++)
//            {
//                h_solution->SConst[i + k*numberMaxConst] = -1;
//            }
//            h_solution->SPAux[k] = 0;
//            h_solution->SMult[0 + k*4] = -1;
//            h_solution->SMult[1 + k*4] = -1;
//            h_solution->SMult[2 + k*4] = -1;
//            h_solution->SMult[3 + k*4] = -1;
//        }
//        free(Coef);
//        free(Coef2);
//        return ;
//
//    }
    for(k=0; k<nRuns; k++)
    {



        violation = 0;
        n1_best = -1;
        n2_best = -1;
        d1_best = -1;
        d2_best = -1;
        for(i=0; i<20; i++)
        {
            Denominator[i]= rand()%maxDenominator + 2;
            Numerator[i] = rand()%(Denominator[i]-1);
        }
        for(mult_1=0; mult_1<20; mult_1++)
        {
            memset(Coef,0,sizeof(int)*h_cut->numberVariables);
            rhs1 = 0;
            for(rest_a = 0; rest_a< numberMaxConst; rest_a++)
            {
                for(i=h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_a] ]; i<h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_a] + 1]; i++)
                {

                    el = h_cut->Elements[i];
                    Coef[el] += h_cut->Coefficients[i] * Numerator[mult_1];
                }
                rhs1 += h_cut->rightSide[ setConstraint[k*numberMaxConst+rest_a] ] * Numerator[mult_1];
                for(mult_2 = 0; mult_2<20; mult_2++)
                {
                    memset(Coef2,0,sizeof(int)*h_cut->numberVariables);
                    value_tes = 0;
                    rhs2 = 0;
                    for(rest_b = rest_a + 1; rest_b<numberMaxConst; rest_b++)
                    {
                        for(j=h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_b] ]; j<h_cut->ElementsConstraints[ setConstraint[k*numberMaxConst + rest_b] + 1]; j++)
                        {
                            el = h_cut->Elements[j];
                            Coef2[el] += h_cut->Coefficients[j] * Numerator[mult_2];
                        }
                        rhs2 += h_cut->rightSide[ setConstraint[k*numberMaxConst + rest_b] ]* Numerator[mult_2];
                    }
                    for(j=0; j<h_cut->numberVariables; j++)
                    {
                        if( (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0 )
                        {
                            aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;

                        }
                        else
                        {
                            aux = (Coef[j]*Denominator[mult_2] + Coef2[j]*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                        }
                        value_tes += aux*h_cut->xAsterisc[j];


//                        aux = Coef[j]<0 ? Coef[j]/Denominator[mult_1] - 1 : Coef[j]/Denominator[mult_1];
//                        value_tes += aux*h_cut->xAsterisc[j];
//                        aux = Coef2[j]<0 ? Coef2[j]/Denominator[mult_2] - 1 : Coef2[j]/Denominator[mult_2];
//                        value_tes += aux*h_cut->xAsterisc[j];
                    }

                    if( (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) < 0)
                    {
                        aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]) - 1;
                    }
                    else
                    {
                        aux = (rhs1*Denominator[mult_2] + rhs2*Denominator[mult_1])/(Denominator[mult_1]*Denominator[mult_2]);
                    }

//                    aux = rhs1<0 ? rhs1/Denominator[mult_1]-1 : rhs1/Denominator[mult_1];
//                    aux +=  rhs2<0 ? rhs2/Denominator[mult_2]-1 : rhs2/Denominator[mult_2];


                    if((value_tes>aux*precision)&&(value_tes-(aux*precision)>violation))
                    {
                        violation = value_tes-(aux*precision);
                        //printf("violation %d\n",violation);
                        n1_best = Numerator[mult_1];
                        d1_best = Denominator[mult_1];
                        n2_best = Numerator[mult_2];
                        d2_best = Denominator[mult_2];
                        qnt_1 = rest_a;
                    }


                }
            }

        }

        if(violation>0)
        {
            for(i=0; i<numberMaxConst; i++)
            {
                h_solution->SConst[i + k*numberMaxConst] = setConstraint[k*numberMaxConst + i];//CPU ja vai ter
            }

            h_solution->SPAux[k] = qnt_1;
            h_solution->SMult[0 + k*4] = n1_best;
            h_solution->SMult[1 + k*4] = d1_best;
            h_solution->SMult[2 + k*4] = n2_best;
            h_solution->SMult[3 + k*4] = d2_best;

        }
        else
        {
            for(i=0; i<numberMaxConst; i++)
            {
                h_solution->SConst[i + k*numberMaxConst] = -1;
            }
            h_solution->SPAux[k] = 0;
            h_solution->SMult[0 + k*4] = -1;
            h_solution->SMult[1 + k*4] = -1;
            h_solution->SMult[2 + k*4] = -1;
            h_solution->SMult[3 + k*4] = -1;
        }


    }
    free(Coef);
    free(Coef2);
}


Cut_gpu* second_phase_runCPU(Cut_gpu *h_cut, Cut_gpu_aux *cut_aux, int numberMaxConst, int nRuns, int maxDenominator, int precision, int szR, double timeLeft)
{

//    printf("FASE 2 in CPU\n");
    int *consR1;
    int *consNR1;
    int *nElemR1;

    double startT = omp_get_wtime();
    double _time = 0;

    Cut_gpu* out_cut_gpu;

    int n_r = 0, n_nr = 0, i,j;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        return h_cut;
    }
    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==RES_R1)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            n_r++;
        }
        else
        {
            n_nr++;
        }
    }
    consR1 = (int*)malloc(sizeof(int)*n_r);
    nElemR1 = (int*)malloc(sizeof(int)*n_r);
    consNR1 = (int*)malloc(sizeof(int)*n_nr);

    n_r = 0;
    n_nr = 0;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consR1);
        free(nElemR1);
        free(consNR1);
        return h_cut;

    }

    for(i=0; i<h_cut->numberConstrains; i++)
    {
        if((h_cut->typeConstraints[i]==RES_RR)||(h_cut->typeConstraints[i]==RES_R1)||(h_cut->typeConstraints[i]==LPC_CGGPU))
        {
            consR1[n_r] = i;
            nElemR1[n_r] = h_cut->ElementsConstraints[i+1] - h_cut->ElementsConstraints[i];
            n_r++;
        }
        else
        {
            if(h_cut->typeConstraints[i]!=LPC_CGGPUR2)
            {
                consNR1[n_nr]=i;
                n_nr++;
            }
        }
    }
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        return h_cut;

    }

    bubble_sort(nElemR1,consR1,n_r);
    int *Similar = returnOrdConstrainsNR(h_cut);
    float *folga = returnFolga(h_cut);
    solutionGpu *h_solution_r2 = allocationStructSolution2(h_cut,numberMaxConst,nRuns);
    int pos_R1 = 0;
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        return h_cut;

    }

    int *setConstraint = (int*)malloc(sizeof(int)*numberMaxConst*nRuns);
    calcSetConstraint(setConstraint, &pos_R1,numberMaxConst, h_cut->numberConstrains, consR1, consNR1, n_r, n_nr, Similar, folga,  nRuns, szR);
    shuffle_Set(setConstraint, numberMaxConst, numberMaxConst*nRuns);
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {
        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        free(setConstraint);
        return h_cut;

    }



    runCPUR2(h_cut, h_solution_r2, numberMaxConst,setConstraint,precision,maxDenominator,nRuns,_time);
    int cont=0;
    //printf("Number constraints: %d\n", h_cut->numberConstrains);
    _time = ( (double) timeLeft - (omp_get_wtime()-startT) );
    if(_time<1)
    {

        free(consNR1);
        free(nElemR1);
        free(consR1);
        free(Similar);
        free(folga);
        free(h_solution_r2);
        free(setConstraint);
        return h_cut;
    }
    for(i=0; i<nRuns; i++)
    {
        if(h_solution_r2->SConst[0 + i*numberMaxConst]!=-1)
        {
            cont++;
        }
    }
    if(cont>0)
    {
         printf("Number of Cuts in the second phase in CPU:%d \n",cont);
        int nT = nRuns;
        int nB = 1;
        //out_cut_gpu = createCutsOfPhaseTwo(h_cut,cut_aux,h_solution_r2,numberMaxConst,cont,precision,nRuns,nT,nB);
        //int nT = nRuns;
        //int nB = 1;
        out_cut_gpu = createCutsStrongPhaseTwo(h_cut,h_solution_r2,numberMaxConst,cont,precision,nRuns,nT,nB);
    }
    else
    {
        free(consR1);
        free(consNR1);
        free(Similar);
        free(folga);
        free(nElemR1);
        free(setConstraint);
        free(h_solution_r2);

        return h_cut;
    }

    free(consR1);
    free(consNR1);
    free(Similar);
    free(folga);
    free(nElemR1);
    free(setConstraint);
    free(h_solution_r2);
    free(h_cut);
    return out_cut_gpu;
}
