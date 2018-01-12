/* at least one solver need to be selected to compile this file */


#ifndef CBC
#ifndef GLPK
#ifndef CPX
#ifndef GRB
#error Compile selecting one solver, i.e. include in compiler parameters: -DCBC or -DGLPK or -DCPX or -DGRB
#endif
#endif
#endif
#endif

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <set>
#include <algorithm>
#include <omp.h>
#include <ctime>
#include <cctype>
#include <string>
#include <vector>
#include <map>
extern "C" {
#include "lp.h"
#ifdef GRB
#include "gurobi_c.h"
#endif
}

#define FILE_NAME_SIZE 1024

#define ERROR( msg ) \
    fprintf( stderr, msg ); \
abort(); \
exit(1);

#ifdef DEBUG
#define LP_CHECK_COL_INDEX( lp, col )\
    if ((col<0)||(col>=lp_cols(lp))) {\
        fprintf( stderr, "ERROR: Invalid column index: %d\n", col );\
        fprintf( stderr, "\tat %s:%d\n", __FILE__, __LINE__ );\
        abort();\
        exit(EXIT_FAILURE);}
#else
#define LP_CHECK_COL_INDEX( lp, col ) {}
#endif

#ifdef DEBUG
#define LP_CHECK_ROW_INDEX( lp, _row )\
    if ((_row<0)||(_row>=lp_rows(lp))) {\
        fprintf( stderr, "ERROR: Invalid row index: %d\n", _row );\
        fprintf( stderr, "\tat %s:%d\n", __FILE__, __LINE__ );\
        abort();\
        exit(EXIT_FAILURE);}
#else
#define LP_CHECK_ROW_INDEX( lp, row ) {}
#endif


#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

static inline bool dbl_equal( const double v1, const double v2 );

static inline bool is_integer( const double v );

using namespace std;

// solver dependent includes and definitions
#ifdef GLPK
#include <glpk.h>
#endif
#ifdef CBC
#ifdef _MSC_VER
#include <OsiSolverInterface.hpp>
#include <OsiCbcSolverInterface.hpp>
#include <OsiClpSolverInterface.hpp>
#include <CoinBuild.hpp>
#include <CglPreProcess.hpp>
#else
#include <coin/OsiSolverInterface.hpp>
#include <coin/OsiClpSolverInterface.hpp>
#include <coin/OsiCbcSolverInterface.hpp>
#include <coin/CoinBuild.hpp>
#include <coin/CglPreProcess.hpp>
#endif
#endif
#ifdef CPX
#include <cplex.h>
static CPXENVptr LPcpxDefaultEnv = NULL;
#endif
#ifdef GRB
static GRBenv *LPgrbDefaultEnv = NULL;
#endif

static bool LPstoreNames = true;

const double EPS=1e-5;

#define INT_NOT_SET INT_MIN

struct _LinearProgram {
    int mipEmphasis;

    char optAsContinuous;
    char mipPreProcess;
    char allInteger; /* if all variables need to be set as integer before solving */

    double obj;
    double bestBound;
    int status;

    int nOptimizations;

    vector< double > *_x;
    vector< double > *_pi;
    vector< double > *_slack;
    vector< double > *_rc;
    vector< double > *_obj;
    vector< int > *_idx;
    vector< double > *_coef;
    vector< int > *_orig;   /* original columns */

    vector< vector< double > > *_savedSol;
    vector< double > *_savedObj;

    map< string, int > *colNameIdx;

    std::vector< std::pair< std::string, double > > *_mipStart;
    std::vector< int > *msIdx;
    std::vector< double > *msVal;

    char silent;

    double solutionTime;

    /* parameters */
    int heurFPPasses;
    int heurProx;
    int maxSeconds;
    int maxSolutions; // exit after found n solutions
    int maxSavedSols;
    int maxNodes;
    int cuts;
    int printMessages;
    double absMIPGap;
    double relMIPGap;
    int parallel;

    /* callback function */
    lp_cb callback_;
    void *data_;

    char solOutFN[256];
    char solInFN[256];

    // defining solver dependent
    // LP storing type
#ifdef GLPK
    glp_prob *_lp;
#endif
#ifdef CBC
    OsiClpSolverInterface *_lp;
    OsiSolverInterface *osiLP;

    /* is a wrapper to an existing osisolverinterface, i.e. should not delete
     * object */
    char justWrapOsi;

    // if model is pre processed then this should be stored
    CglPreProcess *cglPP;

    OsiCuts *cutPool;
    char ownCutPool;
#endif
#ifdef CPX
    CPXLPptr cpxLP;
#endif // CPX
#ifdef GRB
   GRBmodel* lp;
   int tmpRows; // rows not flushed yet
   int tmpCols; // rows not flushed yet
   int nModelChanges;
#endif
};

bool lp_str_is_num( const char *str );

#ifdef GLPK
void lp_config_glpk_params(LinearProgram *lp, glp_iocp *iocp);
#endif
#ifdef CBC
void lp_config_cbc_params(LinearProgram *lp, vector<string> &cbcOpt);
#endif
#ifdef CPX
void lp_check_for_cpx_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine);
void lp_config_cpx_params(LinearProgram *lp);
#endif
#ifdef GRB
void lp_config_grb_params( LinearProgram *lp );
void lp_check_for_grb_error( GRBenv* env, int errorCode, const char *sourceFile, int sourceLine);
#endif

void lp_printRootRelaxInfo(LinearProgram *lp);


#ifdef CBC
/* makes a lp use a predefined OsiSolverInterface (already populated) */
LinearProgram *lp_wrap_osisolver( OsiSolverInterface *osiLP );

// cut generator to accept callbacks in CBC
//
class CglCallback : public CglCutGenerator
{
public:
   CglCallback();

   lp_cb callback_;
   LinearProgram *lp;
   void *data;
   CbcModel *model;

   /// Copy constructor
   CglCallback(const CglCallback& rhs);

   /// Clone
   virtual CglCutGenerator * clone() const;

   virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
         const CglTreeInfo info = CglTreeInfo() );

   virtual ~CglCallback();
private:
};


CglCallback::CglCallback()
    : callback_(NULL),
      lp(NULL),
      data(NULL),
      model(NULL)
{
}

CglCallback::CglCallback(const CglCallback& rhs)
{
    this->callback_ = rhs.callback_;
    this->lp = rhs.lp;
    this->data = rhs.data;
    this->model = rhs.model;
}

CglCutGenerator* CglCallback::clone() const
{
   CglCallback *cglcb = new CglCallback();
   cglcb->callback_ = this->callback_;
   cglcb->lp = this->lp;
   cglcb->data = this->data;
   cglcb->model = this->model;

   return static_cast<CglCutGenerator*>(cglcb);
}

void CglCallback::generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info )
{
    LinearProgram *lp = lp_wrap_osisolver( (OsiSolverInterface *)&si );

    lp->nOptimizations = 0;

    lp->status = LP_ERROR;
    if (lp->osiLP->isProvenOptimal()) {
        /* getting primal and dual solution */
        lp->_x->resize(lp_cols(lp));
        lp->_rc->resize(lp_cols(lp));
        lp->_pi->resize(lp_rows(lp));
        lp->_slack->resize(lp_rows(lp));
        memcpy(&((*(lp->_x))[0]) , lp->osiLP->getColSolution(), sizeof(double)*lp_cols(lp));
        memcpy(&((*(lp->_rc))[0]) , lp->osiLP->getReducedCost(), sizeof(double)*lp_cols(lp));
        memcpy(&((*(lp->_pi))[0]) , lp->osiLP->getRowPrice(), sizeof(double)*lp_rows(lp));
        for ( int i=0 ; (i<lp_rows(lp)) ; ++i )
        {
            double activity = lp->osiLP->getRowActivity()[i];
            double lower = lp->osiLP->getRowLower()[i];
            double upper = lp->osiLP->getRowUpper()[i];
            (*lp->_slack)[i] = std::min( upper-activity, activity-lower );
        }

        lp->obj = lp->osiLP->getObjValue();
        lp->status = LP_OPTIMAL;
        lp->nOptimizations = 1;
    }

    lp->cutPool = &cs;
    this->callback_( lp, LPCB_CUTS, this->model->originalColumns(), this->lp, this->data );
    lp_free(&lp);
}

CglCallback::~CglCallback()
{

}

#endif

/* from names */
void lp_gen_mipstart( LinearProgram *lp )
{
    const int n = lp->_mipStart->size();
    lp->msIdx->clear(); lp->msVal->clear();
    lp->msIdx->reserve( n ); lp->msVal->reserve( n );
    map< string, int > nameIndex; char name[512];

    for (int i=0 ; (i<lp_cols(lp)) ; i++ )
        nameIndex[lp_col_name(lp,i,name)] = i;

    char warned = 0;
    vector< pair< string, double > >::const_iterator mIt;
    for ( mIt=lp->_mipStart->begin() ; (mIt!=lp->_mipStart->end()) ; ++mIt )
    {
        map<string, int>::iterator it = nameIndex.find( mIt->first );
        if ( it == nameIndex.end() )
        {
            if (!warned)
            {
                printf("MIPStart warning: variable %s not found.\n", mIt->first.c_str() );
                warned = 1;
                break;
            }
        }
        else
        {
            lp->msIdx->push_back( it->second );
            lp->msVal->push_back( mIt->second );
        }


    }

}

void lp_initialize(LinearProgram *lp);

LinearProgramPtr lp_create()
{
    LinearProgram *result = (LinearProgramPtr) malloc(sizeof(LinearProgram));
    assert(result);

    lp_initialize(result);

#ifdef GLPK
    result->_lp = glp_create_prob();
    assert(result->_lp);
#endif
#ifdef CBC
    result->_lp   = new OsiClpSolverInterface();
    result->_lp->messageHandler()->setLogLevel(0);
    result->_lp->getModelPtr()->setPerturbation(50);
    result->osiLP = dynamic_cast<OsiSolverInterface *>(result->_lp);
    result->cglPP = NULL;

    if (LPstoreNames)
        result->osiLP->setIntParam(OsiNameDiscipline, 1);
#endif
#ifdef CPX
    int cpxError = 1;

    if (!LPcpxDefaultEnv) {
        LPcpxDefaultEnv = CPXopenCPLEX(&cpxError);
        if (!LPcpxDefaultEnv) {
            fprintf(stderr, "Error opening CPLEX environment. Quiting.\n");
            abort();
        }
    }

    result->cpxLP = CPXcreateprob(LPcpxDefaultEnv, &cpxError, "mip");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif // CPX
#ifdef GRB
    if ( LPgrbDefaultEnv == NULL )
    {
        int status = GRBloadenv(&LPgrbDefaultEnv, NULL);
        if (status)
        {
            fprintf( stderr, "ERROR: Gurobi environment could not be loaded. Check license.\n");
            abort();
        }
    }

    result->lp = NULL;
    int status = GRBnewmodel( LPgrbDefaultEnv, &result->lp, "", 0, NULL, NULL, NULL, NULL, NULL);
    if (status)
    {
        fprintf( stderr, "ERROR: Could not create Gurobi Model.\n");
        abort();
    }

    result->tmpRows = 0;
    result->tmpCols = 0;
    result->nModelChanges = 0;
#endif

    result->allInteger = 0;

    return result;
}

char getFileType(const char *fileName)
{
    if ( strstr(fileName, ".mps") || strstr(fileName, ".MPS") || strstr(fileName, ".mps.gz") || strstr(fileName, ".MPS.GZ") )
        return 'M';

    return 'L';
}

void lp_fill_col_name_index( LinearProgram *lp )
{
    lp->colNameIdx->clear();
    char colName[512];

    for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
        (*lp->colNameIdx)[lp_col_name(lp, i, colName)] = i;
}

#ifdef CBC
LinearProgram *lp_wrap_osisolver( OsiSolverInterface *osiLP )
{
    LinearProgram *lp = (LinearProgramPtr) calloc(sizeof(LinearProgram), 1);
    assert(lp);

    lp_initialize(lp);

    lp->justWrapOsi = 1;
    lp->osiLP = osiLP;
    lp_fill_col_name_index( lp );

    return lp;
}
#endif

void lp_read(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);

#ifdef CPX
    int cpxError;
    cpxError = CPXreadcopyprob(LPcpxDefaultEnv, lp->cpxLP, fileName, NULL);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    lp_fill_col_name_index( lp );
    return;
#endif
#ifdef GRB
    int grbError = GRBreadmodel(LPgrbDefaultEnv, fileName, &(lp->lp) );
    lp_check_for_grb_error(LPgrbDefaultEnv, grbError, __FILE__, __LINE__);
    lp_fill_col_name_index( lp );
    return;
#endif

    switch (getFileType(fileName)) {
        case 'L':
#ifdef GLPK
            glp_read_lp(lp->_lp, NULL, fileName);
#endif
#ifdef CBC
            lp->osiLP->readLp(fileName);
#endif
            break;
        case 'M':
#ifdef GLPK
            glp_read_mps(lp->_lp, GLP_MPS_FILE, NULL, fileName);
#endif
#ifdef CBC
            lp->osiLP->readMps(fileName);
#endif
            break;
    }
    lp_fill_col_name_index( lp );
}

void lp_write_lp(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);
    assert( fileName != NULL );
    assert( strlen(fileName) );

#ifdef CPX
    int cpxError;
    char fName[FILE_NAME_SIZE];
    strcpy( fName, fileName );
    if ( (!strstr(fName,".lp")) && (!strstr(fName,".LP")) )
        strcat( fName, ".lp" );
    cpxError = CPXwriteprob(LPcpxDefaultEnv, lp->cpxLP, fName, "LP");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return;
#endif
#ifdef GRB
    if (lp->tmpRows>0)
    {
        int grbError;
        grbError = GRBupdatemodel(lp->lp);
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

        lp->tmpRows = 0;
    }

    char fName[256];
    strcpy( fName, fileName );
    if ( strstr(fName, ".lp")==0 )
        strcat( fName, ".lp" );

    int grbError = GRBwrite( lp->lp, fName );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return;
#endif

    char fileType = getFileType(fileName);

    switch (fileType) {
        case 'L':
#ifdef GLPK
            glp_write_lp(lp->_lp, NULL, fileName);
#endif
#ifdef CBC
            {
                char outFile[256];
                strcpy(outFile, fileName);
                char *s = NULL;
                if ((s = strstr(outFile, ".lp"))) {
                    if (s != outFile) // not at the start
                        *s = '\0';
                }
                lp->osiLP->writeLp(outFile);
            }
#endif

            break;
        case 'M':
#ifdef GLPK
            glp_write_mps(lp->_lp, GLP_MPS_FILE, NULL, fileName);
#endif
#ifdef CBC
            lp->osiLP->writeMps(fileName);
#endif

            break;
    }
}

void lp_set_direction(LinearProgram *lp, const char direction)
{
    assert(lp != NULL);

    switch (direction) {
        case LP_MIN:
#ifdef GLPK
            glp_set_obj_dir(lp->_lp, GLP_MIN);
#endif
#ifdef CBC
            lp->osiLP->setObjSense(1.0);
#endif
#ifdef CPX
            CPXchgobjsen(LPcpxDefaultEnv, lp->cpxLP, CPX_MIN);
#endif
#ifdef GRB
            {
                int grbError = GRBsetintattr( lp->lp, GRB_INT_ATTR_MODELSENSE, 1 );
                lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
                grbError = GRBupdatemodel(lp->lp);
                lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
                return;
            }
#endif
            break;
        case LP_MAX:
#ifdef GLPK
            glp_set_obj_dir(lp->_lp, GLP_MAX);
#endif
#ifdef CBC
            lp->osiLP->setObjSense(-1.0);
#endif
#ifdef GRB
            {
                int grbError = GRBsetintattr( lp->lp, GRB_INT_ATTR_MODELSENSE, -1 );
                lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
                grbError = GRBupdatemodel(lp->lp);
                lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
                return;
            }
#endif
#ifdef CPX
            CPXchgobjsen(LPcpxDefaultEnv, lp->cpxLP, CPX_MAX);
#endif
            break;
        default:
            fprintf(stderr, "Unknow optimization direction: %c\n", direction);
            break;
    }
}

int lp_get_direction(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef GLPK
    if (glp_get_obj_dir(lp->_lp) ==  GLP_MAX)
        return LP_MAX;
    return LP_MIN;
#endif
#ifdef CBC
    if ((fabs(lp->osiLP->getObjSense() - 1.0)) <= EPS)
        return LP_MIN;
    return LP_MAX;
#endif
#ifdef GRB
    int dir = INT_MAX;
    int grbError = GRBgetintattr( lp->lp, GRB_INT_ATTR_MODELSENSE, &dir );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    switch (dir)
    {
        case 1:
            return LP_MIN;
        case -1:
            return LP_MAX;
        default:
        {
            fprintf( stderr, "Invalid optimization directtion: %d\n", dir );
            abort();
        }
    }
#endif
#ifdef CPX
    switch (CPXgetobjsen(LPcpxDefaultEnv, lp->cpxLP)) {
        case CPX_MIN:
            return LP_MIN;
            break;
        case CPX_MAX:
            return LP_MAX;
            break;
        default:
            fprintf(stderr, "Invalid value for CPXgetobjsen.\n");
            abort();
    }
#endif
}

#ifdef DEBUG
static void lp_validate_row_data(LinearProgram *lp, const int nz,  int *indexes, double *coefs, const char *name, char sense, const double rhs)
{
    assert( lp != NULL );
    assert( nz>=1 );
    assert( indexes != NULL );
    assert( coefs != NULL );
    assert( nz <= lp_cols(lp) );

    /* checking indexes */
    for (int i = 0 ; (i < nz) ; ++i) {
        LP_CHECK_COL_INDEX(lp, indexes[i]);
        for ( int j=i+1 ; j<nz ; ++j )
        {
            if ( indexes[i]==indexes[j] )
            {
#ifdef GRB
                    if ( lp->tmpCols || lp->tmpRows || lp->nModelChanges )
                    {
                        int grbError = GRBupdatemodel( lp->lp );
                        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
                        lp->tmpCols = lp->tmpRows = lp->nModelChanges = 0;
                    }
#endif
                char varName[256];
                lp_col_name( lp, indexes[i], varName );
                fprintf( stderr, "ERROR adding row in Linear Program: variable (%s) appears more than once in a constraint (%s).\n", varName, name );
                abort();
            }
        } // other nz
    } // nzs
}
#endif

void lp_add_cut( LinearProgram *lp, int nz, int *cutIdx, double *cutCoef, const char *name, char sense, double rhs )
{
#ifdef DEBUG
    lp_validate_row_data( lp, nz, cutIdx, cutCoef, name, sense, rhs );
#endif

    sense = toupper(sense);

#ifdef CBC
    OsiRowCut orCut;
    orCut.setRow( nz, cutIdx, cutCoef, false );

    double rLB = 0.0, rUB = COIN_DBL_MAX;
    switch (sense) {
        case 'E':
            rLB = rhs;
            rUB = rhs;
            break;
        case 'L':
            rLB = -COIN_DBL_MAX;
            rUB = rhs;
            break;
        case 'G':
            rLB = rhs;
            rUB = COIN_DBL_MAX;
            break;
    }
    orCut.setLb( rLB );
    orCut.setUb( rUB );

    lp->cutPool->insert( orCut );
#endif
#ifdef CPX
    lp_add_row( lp, nz, cutIdx, cutCoef, name, sense, rhs );
#endif
#ifdef GRB
    lp_add_row( lp, nz, cutIdx, cutCoef, name, sense, rhs );
#endif
}

void lp_add_row(LinearProgram *lp, const int nz,  int *indexes, double *coefs, const char *name, char sense, const double rhs)
{
#ifdef DEBUG
    lp_validate_row_data( lp, nz, indexes, coefs, name, sense, rhs );

    static std::set< string > rowNames;
    static int nWarnings = 0;
    if ( rowNames.find(name) != rowNames.end() )
    {
        if ( nWarnings++ == 0 )
            fprintf( stdout, "[warning] inserting rowName %s twice.\n", name );
        rowNames.insert( name );
    }

    // checkinf if there are coefficients zero
    for ( int i=0 ; i<nz ; ++i )
        if (fabs(coefs[i])<1e-13)
        {
            char colName[256];
            lp_col_name( lp, indexes[i], colName );
            printf("warning: coefficient %g for variable %s ignored while adding constraint %s.\n", coefs[i], colName, name );
        }
#endif

    sense = toupper(sense);

    char strsense[2]; sprintf(strsense, "%c", sense);
    if (!strstr("GLE", strsense)) {
        fprintf(stderr, "LP: cannot handle sense %c.\n", sense);
        abort();
    }

#ifdef CPX
    int cpxError;

    int matBeg[] = { 0, nz };

    cpxError = CPXaddrows(LPcpxDefaultEnv, lp->cpxLP, 0, 1, nz, &rhs, &sense, matBeg, indexes, coefs, NULL, (char **) &name);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return;
#endif // CPX
#ifdef GRB
    // making sure that gurobi defines are respected:
    switch (sense)
    {
        case 'E':
            sense = GRB_EQUAL;
            break;
        case 'L':
            sense = GRB_LESS_EQUAL;
            break;
        case 'G':
            sense = GRB_GREATER_EQUAL;
            break;
    }

    lp->tmpRows++;

    int grbError = GRBaddconstr( lp->lp, nz, indexes, coefs, sense, rhs, name );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif
#ifdef GLPK
    int currRow = lp_rows(lp) + 1;

    glp_add_rows(lp->_lp, 1);
    glp_set_row_name(lp->_lp, currRow, name);
    switch (sense) {
        case 'L' :
            glp_set_row_bnds(lp->_lp, currRow, GLP_UP, 0.0, rhs);
            break;
        case 'G' :
            glp_set_row_bnds(lp->_lp, currRow, GLP_LO, rhs, 0.0);
            break;
        case 'E' :
            glp_set_row_bnds(lp->_lp, currRow, GLP_FX, rhs, 0.0);
            break;
    }

    register int *endIdx = indexes + nz;
    register int *idxPtr;
    for (idxPtr = indexes ; (idxPtr < endIdx) ; idxPtr++)
        (*idxPtr)++;

    glp_set_mat_row(lp->_lp, currRow, nz, indexes - 1, coefs - 1);

    // restoring indexes
    for (idxPtr = indexes ; (idxPtr < endIdx) ; idxPtr++)
        (*idxPtr)--;
#endif
#ifdef CBC
    int currRow = lp_rows(lp);

    double rLB = 0.0, rUB = COIN_DBL_MAX;
    switch (sense) {
        case 'E':
            rLB = rhs;
            rUB = rhs;
            break;
        case 'L':
            rLB = -COIN_DBL_MAX;
            rUB = rhs;
            break;
        case 'G':
            rLB = rhs;
            rUB = COIN_DBL_MAX;
            break;
    }
    lp->osiLP->addRow(nz, indexes, coefs, rLB, rUB);
    lp->osiLP->setRowName(currRow, name);
#endif
}

void lp_add_cols(LinearProgram *lp, const int count, double *obj, double *lb, double *ub, char *integer, char **name)
{
#ifdef DEBUG
    assert( lp != NULL );
    assert( count >=1  );
    assert( name != NULL );

    if ( (lb) && (ub) )
        for ( int i=0 ; (i<count) ; i++ )
            assert(lb[i] <= ub[i]);
#endif

#ifdef GRB
    vector< char > type( count, GRB_CONTINUOUS );

    double *_lb = lb;
    double *_ub = ub;

    vector< double > rLB;
    vector< double > rUB;
    if (!_lb) {
        rLB.resize(count, 0.0);
        _lb = &rLB[0];
    }
    if (!_ub) {
        rUB.resize(count, GRB_INFINITY );
        _ub = &rUB[0];
    }

    for (int i = 0 ; (i < count) ; i++)
    {
        if (_ub[i] > GRB_INFINITY)
            _ub[i] = GRB_INFINITY;
    }

    if (integer) {
        for (int i = 0 ; (i < count) ; ++i)
            if (integer[i]) {
                _lb[i] = floor(_lb[i] + 0.5);

                if (_ub[i] != GRB_INFINITY)
                    _ub[i] = floor(_ub[i] + 0.5);

                if ((fabs(_lb[i]) < EPS) && (fabs(_ub[i] - 1.0) < EPS)) {
                    type[i] = GRB_BINARY;
                }
                else {
                    type[i] = GRB_INTEGER;
                }
            }
    }

    vector< int > vbeg( count+1 , 0 );
    int vind[2] = { 0, 0 };
    double vval[2] = { 0, 0 };
    int grbError = GRBaddvars( lp->lp, count, 0, &vbeg[0], vind, vval, obj, lb, ub, &type[0], name );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    grbError = GRBupdatemodel(lp->lp);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif
#ifdef CPX
    int cpxError;

    vector< char > type( count, CPX_CONTINUOUS );

    double *_lb = lb;
    double *_ub = ub;

    vector< double > rLB;
    vector< double > rUB;
    if (!_lb) {
        rLB.resize(count, 0.0);
        _lb = &rLB[0];
    }
    if (!_ub) {
        rUB.resize(count, CPX_INFBOUND);
        _ub = &rUB[0];
    }

    for (int i = 0 ; (i < count) ; i++) {
        if (_ub[i] > CPX_INFBOUND)
            _ub[i] = CPX_INFBOUND;
    }


    if (integer) {
        for (int i = 0 ; (i < count) ; ++i)
            if (integer[i]) {
                _lb[i] = floor(_lb[i] + 0.5);

                if (_ub[i] != CPX_INFBOUND)
                    _ub[i] = floor(_ub[i] + 0.5);

                if ((fabs(_lb[i]) < EPS) && (fabs(_ub[i] - 1.0) < EPS)) {
                    type[i] = CPX_BINARY;
                }
                else {
                    type[i] = CPX_INTEGER;
                }
            }
    }

    cpxError = CPXnewcols(LPcpxDefaultEnv, lp->cpxLP, count, obj, _lb, _ub, &type[0] , name);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif // CPX
#ifdef GLPK
    register int j, cols, currCol;
    cols = lp_cols(lp);

    glp_add_cols(lp->_lp, count);

    if (name)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            glp_set_col_name(lp->_lp, currCol, name[j]);

    if (obj)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            glp_set_obj_coef(lp->_lp, currCol, obj[j]);

    if (integer)
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            if (integer[j]) {
                if ((fabs(ub[j] - 1.0) <= EPS))
                    glp_set_col_kind(lp->_lp, currCol, GLP_BV);
                else
                    glp_set_col_kind(lp->_lp, currCol, GLP_IV);
            }

    for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++) {
        if (((ub) && (lb))) {
            if ((lb[j] != -DBL_MAX) && (ub[j] != DBL_MAX))
            {
                if ( fabs(lb[j]-ub[j]) <= 1e-10 )
                    glp_set_col_bnds(lp->_lp, currCol, GLP_FX, lb[j] , ub[j]);
                else
                    glp_set_col_bnds(lp->_lp, currCol, GLP_DB, lb[j] , ub[j]);
            }
            else if ((ub[j] == DBL_MAX) && (lb[j] != -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, lb[j] , ub[j]);
            else if ((ub[j] != DBL_MAX) && (lb[j] == -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_UP, lb[j] , ub[j]);
            else if ((ub[j] == DBL_MAX) && (lb[j] == -DBL_MAX))
                glp_set_col_bnds(lp->_lp, currCol, GLP_FR, -DBL_MAX , DBL_MAX);
        }   // both LB and UB informed
        else {
            if ((!lb) && (!ub)) {
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
            }  // no LB and UB informed
            else {
                if (ub) {
                    if (ub[j] == DBL_MAX)
                        glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
                    else
                        glp_set_col_bnds(lp->_lp, currCol, GLP_DB, 0.0 , ub[j]);
                }
                else if (lb) {
                    if (lb[j] == -DBL_MAX)
                        glp_set_col_bnds(lp->_lp, currCol, GLP_FR, -DBL_MAX , DBL_MAX);
                    else
                        glp_set_col_bnds(lp->_lp, currCol, GLP_DB, lb[j] , DBL_MAX);
                }
            }  // only LB or UB informed
        }
    } // all columns

    // configuring positive, continuous variables
    // (for binary vars it is not necessary to inform
    // bounds)
    if ((integer) && (!ub))
        for (currCol = cols + 1, j = 0 ; (j < count) ; j++, currCol++)
            if (!integer[j])
                glp_set_col_bnds(lp->_lp, currCol, GLP_LO, 0.0 , DBL_MAX);
#endif
#ifdef CBC
    int cols = lp->osiLP->getNumCols();

    {
        vector< int > starts(count + 1, 0);
        lp->osiLP->addCols(count, &starts[0], NULL, NULL, lb, ub, obj);
    }


    if (integer) {
        for (int j = 0 ; (j < count) ; j++)
            if (integer[j])
                lp->osiLP->setInteger(cols + j);
    }

    for (int j = 0 ; (j < count) ; j++)
        lp->osiLP->setColName(cols + j, name[j]);
#endif
    //  updating column name index
    {
        int idxCol = lp_cols(lp) - count;
        for (int j = 0 ; (j < count) ; j++)
        {
            map< string, int >::const_iterator mIt;
            mIt = (*lp->colNameIdx).find(name[j]);
            if ( mIt != (*lp->colNameIdx).end() )
            {
                fprintf( stderr, "ERROR: trying to add variable with repeated name: %s\n", name[j]);
                abort();
            }

            (*lp->colNameIdx)[name[j]] = idxCol++;
        }
    }
}

void lp_add_cols_same_bound(LinearProgram *lp, const int count, double *obj, double lb, double ub, char *integer, char **name)
{
    assert(lp != NULL);

    vector< double > vlb(count, lb);
    vector< double > vub(count, ub);
    lp_add_cols(lp, count, obj, &vlb[0], &vub[0], integer, name);
}

void lp_add_bin_cols( LinearProgram *lp, const int count, double *obj, char **name )
{
    vector< double > vlb(count, 0.0 );
    vector< double > vub(count, 1.0 );
    vector< char > integer( count, 1 );
    lp_add_cols( lp, count, obj, &vlb[0], &vub[0], &(integer[0]), name );
}

const double *lp_obj_coef( LinearProgram *lp )
{
    assert(lp);

#ifdef CBC
    OsiClpSolverInterface *osilp = lp->_lp;

    return osilp->getObjCoefficients();
#endif
#ifdef GLPK
    {
        int i;
        glp_prob *glp = lp->_lp;
        lp->_obj->resize( lp_cols(lp) );
        for ( i=0 ; (i<lp_cols(lp)) ; i++ )
            (*lp->_obj)[i] = glp_get_obj_coef( glp, i+1 );
    }
    return &((*lp->_obj)[0]);
#endif
#ifdef CPX
    {
        lp->_obj->resize( lp_cols(lp) );
        int cpxError =  CPXgetobj( LPcpxDefaultEnv, lp->cpxLP, &((*lp->_obj)[0]), 0, lp_cols(lp)-1 );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
    return &((*lp->_obj)[0]);
#endif
#ifdef GRB
    lp->_obj->resize( lp_cols(lp) );
    double *obj = &((*lp->_obj)[0]);
    int grbError = GRBgetdblattrarray( lp->lp, GRB_DBL_ATTR_OBJ, 0, lp_cols(lp) , obj );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return obj;
#endif

    return NULL;
}

int lp_cols(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef GLPK
    return glp_get_num_cols(lp->_lp);
#endif
#ifdef CBC
    return lp->osiLP->getNumCols();
#endif
#ifdef GRB
    int numCols = INT_MAX;
    int grbError = GRBgetintattr(lp->lp, GRB_INT_ATTR_NUMVARS, &numCols);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return numCols;
#endif
#ifdef CPX
    return CPXgetnumcols(LPcpxDefaultEnv, lp->cpxLP);
#endif
}

int lp_rows(LinearProgram *lp)
{
    assert( lp != NULL );

#ifdef GLPK
    return glp_get_num_rows(lp->_lp);
#endif
#ifdef CBC
    return lp->osiLP->getNumRows();
#endif
#ifdef GRB
    int numRows = INT_MAX;
    int grbError = GRBgetintattr(lp->lp, GRB_INT_ATTR_NUMCONSTRS, &numRows);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return numRows + lp->tmpRows;
#endif
#ifdef CPX
    return CPXgetnumrows(LPcpxDefaultEnv, lp->cpxLP);
#endif
}


char lp_is_integer(LinearProgram *lp, const int j)
{
    assert( lp != NULL );
    LP_CHECK_COL_INDEX(lp, j);

#ifdef GLPK
    switch (glp_get_col_kind(lp->_lp, j + 1)) {
        case GLP_IV:
        case GLP_BV:
            return 1;
            break;
    }

    return 0;
#endif
#ifdef CBC
    return lp->osiLP->isInteger(j);
#endif
#ifdef CPX
    if (CPXgetprobtype(LPcpxDefaultEnv, lp->cpxLP) == CPXPROB_LP)
        return 0;

    char colType[2];
    int cpxError = CPXgetctype(LPcpxDefaultEnv, lp->cpxLP, colType, j, j);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return ((colType[0] == CPX_BINARY) || (colType[0] == CPX_INTEGER));
#endif
#ifdef GRB
    char vType = 0;
    int grbError = GRBgetcharattrelement( lp->lp, GRB_CHAR_ATTR_VTYPE, j, &vType );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    switch (vType)
    {
        case GRB_INTEGER :
            return 1;
        case GRB_BINARY :
            return 1;
        default :
            return 0;
    }

    return 0;
#endif
}

char lp_isMIP(LinearProgram *lp)
{
    int nCols = lp_cols(lp);
    int j;

    for (j = 0 ; (j < nCols) ; j++)
        if (lp_is_integer(lp, j))
            return 1;

    return 0;
}

int lp_optimize_as_continuous(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->optAsContinuous = 1;
    int status = lp_optimize(lp);
    lp->optAsContinuous = 0;

    return status;
}

int lp_optimize(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->solutionTime = 0.0;

    time_t startT; time(&startT);

#ifdef GRB
    if ( (lp->tmpRows>0) || (lp->tmpCols>0) || (lp->nModelChanges) )
    {
        int grbError;
        grbError = GRBupdatemodel(lp->lp);
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

        lp->tmpRows = 0;
        lp->tmpCols = 0;
        lp->nModelChanges = 0;
    }
#endif

    int isMIP = 0;

    if (!lp->optAsContinuous)
        isMIP = lp_isMIP(lp);

    lp->_x->resize(lp_cols(lp));
    lp->_rc->resize(lp_cols(lp));
    lp->_obj->resize(lp_cols(lp));
    lp->_pi->resize(lp_rows(lp));
    lp->_slack->resize(lp_rows(lp));
    lp->_savedSol->clear();
    lp->_savedObj->clear();
    lp->obj = DBL_MAX;
    lp->bestBound = DBL_MAX;
    lp->status = LP_ERROR;

    fill(lp->_x->begin(), lp->_x->end(), DBL_MAX);
    fill(lp->_rc->begin(), lp->_rc->end(), DBL_MAX);
    fill(lp->_obj->begin(), lp->_obj->end(), DBL_MAX);
    fill(lp->_pi->begin(), lp->_pi->end(), DBL_MAX);
    fill(lp->_slack->begin(), lp->_slack->end(), DBL_MAX);

    lp->nOptimizations++;

    /* error handling */
    char errorMsg[256] = "";
    int errorLine = -1;

    /* check if problem will be transformed in all integer */
    if (lp->allInteger)
    {
        int cInt, cBin, cCont;
        lp_cols_by_type( lp, &cBin, &cInt, &cCont );
        if (cBin+cInt < lp_cols(lp))
        {
            vector< int > idx;
            for ( int i=0 ; (i<lp_cols(lp)) ; ++i )
                if (!lp_is_integer(lp,i))
                    idx.push_back(i);

            lp_set_integer( lp, idx.size(), &idx[0] );
        }
    } /* transforming into pure IP */

#ifdef GLPK
    {
        // solving linear relaxation
        int status = 1;

        glp_smcp parm; glp_init_smcp(&parm);

        if (lp->silent)
            parm.msg_lev = GLP_MSG_OFF;

        status = glp_simplex(lp->_lp, &parm);

        // checking error type
        if (status) {
            switch (status) {
                case GLP_EBADB :
                    sprintf(errorMsg,
                            "Unable to start the search, because the initial basis specified in the problem object is\
                            invalid—the number of basic (auxiliary and structural) variables is not the same as the\
                            number of rows in the problem object.");
                    errorLine = __LINE__;
                    break;
                case GLP_ESING :
                    sprintf(errorMsg,
                            "Unable to start the search, because the basis matrix corresponding\
                            to the initial basis is singular within the working precision.");
                    errorLine = __LINE__;
                    break;
                case GLP_ECOND :
                    sprintf(errorMsg,
                            "Unable to start the search, because the basis matrix corresponding to the initial basis is\
                            ill-conditioned, i.e. its condition number is too large.");
                    errorLine = __LINE__;
                    break;
                case GLP_EBOUND :
                    sprintf(errorMsg,
                            "Unable to start the search, because some double-bounded\
                            (auxiliary or structural) variables have incorrect bounds.");
                    errorLine = __LINE__;
                    break;
                case GLP_EFAIL :
                    sprintf(errorMsg,
                            "The search was prematurely terminated due to the solver failure.");
                    errorLine = __LINE__;
                    break;
                case GLP_EOBJLL :
                    sprintf(errorMsg,
                            "The search was prematurely terminated, because the objective function being maximized has reached its\
                            lower limit and continues decreasing (the dual simplex only).");
                    errorLine = __LINE__;
                    break;
                case GLP_EOBJUL :
                    sprintf(errorMsg,
                            "The search was prematurely terminated, because the objective function being minimized has reached its\
                            upper limit and continues increasing (the dual simplex only).");
                    errorLine = __LINE__;
                    break;
                case GLP_EITLIM :
                    sprintf(errorMsg,
                            "The search was prematurely terminated, because the simplex iteration limit has been exceeded.");
                    errorLine = __LINE__;
                    break;
                case GLP_ETMLIM :
                    sprintf(errorMsg,
                            "The search was prematurely terminated, because the time limit has been exceeded.");
                    errorLine = __LINE__;
                    break;
                case GLP_ENOPFS :
                    sprintf(errorMsg,
                            "The LP problem instance has no primal feasible solution (only if the LP presolver is used).");
                    errorLine = __LINE__;
                    break;
                case GLP_ENODFS :
                    sprintf(errorMsg,
                            "The LP problem instance has no dual feasible solution (only if the LP presolver is used).");
                    errorLine = __LINE__;
                    break;
            }

            goto OPTIMIZATION_ERROR;
        }
        else {
            status = glp_get_status(lp->_lp);

            switch (status) {
                case GLP_OPT:
                    if (!isMIP) {
                        for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
                        {
                            lp->_x->at(i) = glp_get_col_prim(lp->_lp, i + 1);
                            lp->_rc->at(i) = glp_get_col_dual(lp->_lp, i + 1);
                        }
                        for (int i = 0 ; (i < lp_rows(lp)) ; ++i)
                        {
                            lp->_pi->at(i) = glp_get_row_dual(lp->_lp, i + 1);
                            double activity = glp_get_row_prim( lp->_lp, i + 1 );
                            double slack = DBL_MAX;
                            switch (glp_get_row_type(lp->_lp, i + 1))
                            {
                                case GLP_LO:
                                {
                                    slack = activity - glp_get_row_lb(lp->_lp, i + 1);
                                    break;
                                }
                                case GLP_UP:
                                {
                                    slack = glp_get_row_ub(lp->_lp, i + 1) - activity;
                                    break;
                                }
                                case GLP_DB:
                                {
                                    double s1 = activity - glp_get_row_lb(lp->_lp, i + 1);
                                    double s2 = glp_get_row_ub(lp->_lp, i + 1) - activity;
                                    slack = std::min( s1, s2 );
                                    break;
                                }
                                case GLP_FX:
                                {
                                    slack = 0.0;
                                    break;
                                }
                            }
                            lp->_slack->at(i) = slack;
                        }

                        lp->obj = glp_get_obj_val(lp->_lp);
                        lp->status = LP_OPTIMAL;
                        goto OPTIMIZATION_CONCLUDED;
                    }
                    break;
                case GLP_INFEAS:
                case GLP_NOFEAS:
                    lp->status = LP_INFEASIBLE;
                    goto OPTIMIZATION_CONCLUDED;
                case GLP_UNBND:
                    goto OPTIMIZATION_CONCLUDED;
                    lp->status = LP_UNBOUNDED;
                default:
                    sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_simplex at file %s\n", __FILE__);
                    errorLine = __LINE__;
                    goto OPTIMIZATION_ERROR;
            }
        }

        if (isMIP) {
            {
                glp_iocp ioParams;
                lp_config_glpk_params(lp, &ioParams);
                status = glp_intopt(lp->_lp, &ioParams);
            }

            switch (status) {
                case GLP_EFAIL :
                case GLP_EROOT :
                case GLP_EBOUND :
                    sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_intopt at file %s\n", __FILE__);
                    errorLine = __LINE__;
                    goto OPTIMIZATION_ERROR;
            }

            switch (glp_mip_status(lp->_lp)) {
                case GLP_OPT:
                    lp->status = LP_OPTIMAL;
                    goto GLPK_GET_MIP_SOLUTION;
                case GLP_FEAS:
                    lp->status = LP_FEASIBLE;
                    goto GLPK_GET_MIP_SOLUTION;
                case GLP_NOFEAS:
                    lp->status = LP_INFEASIBLE;
                    goto OPTIMIZATION_CONCLUDED;
                case GLP_UNDEF:
                    lp->status = LP_NO_SOL_FOUND;
                    goto OPTIMIZATION_CONCLUDED;
                default:
                    sprintf(errorMsg, "\n\nGLPK ERROR CALLING glp_mip_status at file %s\n", __FILE__);
                    errorLine = __LINE__;
                    goto OPTIMIZATION_ERROR;
            }

GLPK_GET_MIP_SOLUTION:
            for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
                lp->_x->at(i) = glp_mip_col_val(lp->_lp, i + 1);

            lp->obj = glp_mip_obj_val(lp->_lp);
            goto OPTIMIZATION_CONCLUDED;
        }
    }
#endif
#ifdef CBC
    bool deleteLP = false;
    OsiSolverInterface *linearProgram = NULL;

    {
        lp->status = LP_ERROR;

        if ( lp->nOptimizations == 1 )
            lp->_lp->getModelPtr()->setPerturbation(50);

        if (lp_isMIP(lp)) {
            linearProgram = lp->osiLP->clone();
            deleteLP = true;
        }
        else
            linearProgram = lp->osiLP;


        if ( lp->maxSeconds != INT_NOT_SET )
        {
            ClpSimplex *clp = lp->_lp->getModelPtr();
            clp->setMaximumSeconds( lp->maxSeconds );
        }

        if (lp_isMIP(lp))
            linearProgram->initialSolve();
        else {
            if (lp->nOptimizations >= 2)
                linearProgram->resolve();
            else
                linearProgram->initialSolve();
        }

        if (linearProgram->isAbandoned()) {
            sprintf(errorMsg, "Numerical difficulties, linear program abandoned.\n");
            errorLine = __LINE__;
            goto CBC_OPTIMIZATION_ERROR;
        }
        if ((linearProgram->isProvenPrimalInfeasible()) || (linearProgram->isProvenDualInfeasible())) {
            lp->status = LP_INFEASIBLE;
            goto CBC_OPTIMIZATION_CONCLUDED;
        }
        if (linearProgram->isIterationLimitReached()) {
            sprintf(errorMsg, "Iteration limit for linear program reached, abandoning.\n");
            errorLine = __LINE__;
            goto CBC_OPTIMIZATION_ERROR;
        }
        if (linearProgram->isPrimalObjectiveLimitReached()) {
            sprintf(errorMsg, "Primal objective limit reached, abandoning.\n");
            errorLine = __LINE__;
            goto CBC_OPTIMIZATION_ERROR;
        }
        if (linearProgram->isDualObjectiveLimitReached()) {
            sprintf(errorMsg, "Dual objective limit reached, abandoning.\n");
            errorLine = __LINE__;
            goto CBC_OPTIMIZATION_ERROR;
        }

        if ((!isMIP) && (linearProgram->isProvenOptimal())) {
            memcpy(&((*(lp->_x))[0]) , linearProgram->getColSolution(), sizeof(double)*lp_cols(lp));
            memcpy(&((*(lp->_pi))[0]), linearProgram->getRowPrice(), sizeof(double)*lp_rows(lp));
            memcpy(&((*(lp->_rc))[0]), linearProgram->getReducedCost(), sizeof(double)*lp_cols(lp));
            /* computing slack */
            for ( int i=0 ; (i<lp_rows(lp)) ; ++i )
            {
                double activity = lp->osiLP->getRowActivity()[i];
                double lower = lp->osiLP->getRowLower()[i];
                double upper = lp->osiLP->getRowUpper()[i];

                (*lp->_slack)[i] = std::min( upper-activity, activity-lower );
            }

            lp->obj = linearProgram->getObjValue();

            lp->status = LP_OPTIMAL;
            goto CBC_OPTIMIZATION_CONCLUDED;
        }

        if (isMIP) {
            // Pass to Cbc initialize defaults
            CbcModel modelA(*linearProgram);
            if (lp->_mipStart->size())
                modelA.setMIPStart( *(lp->_mipStart) );

            CglCallback *cglCB = NULL;
            if ( lp->callback_ != NULL )
            {
                cglCB = new CglCallback();
                cglCB->callback_ = lp->callback_;
                cglCB->lp = lp;
                cglCB->data = lp->data_;
                cglCB->model = &modelA;
                modelA.addCutGenerator( cglCB, -1, "Callback" );
            }

            CbcModel *model = &modelA;
            {
                // filling options and calling solver
#define STR_OPT_SIZE 256
                vector<string> cbcP;
                char **cbcOptStr = NULL;
                CbcMain0(modelA);
                lp_config_cbc_params(lp, cbcP);
                cbcOptStr = (char **) malloc(sizeof(char *)*cbcP.size()); assert(cbcOptStr);
                cbcOptStr[0] = (char *)malloc(sizeof(char) * cbcP.size() * STR_OPT_SIZE); assert(cbcOptStr[0]);
                for (int i = 1 ; (i < (int)cbcP.size()) ; i++) cbcOptStr[i] = cbcOptStr[i - 1] + STR_OPT_SIZE;
                for (int i = 0 ; (i < (int)cbcP.size()) ; i++) strncpy(cbcOptStr[i], cbcP[i].c_str(), STR_OPT_SIZE);

//                printf("calling CBC with options: %s %d\n", (const char **)cbcOptStr, cbcP.size() );

                CbcMain1(cbcP.size(), (const char **)cbcOptStr, modelA);

                if (cglCB)
                    delete cglCB;
#undef STR_OPT_SIZE
                free(cbcOptStr[0]); free(cbcOptStr); cbcOptStr = NULL;
            }

            lp->status = LP_NO_SOL_FOUND;

            if (model->isAbandoned()) {
                sprintf(errorMsg, "Model isAbandoned()\n");
                errorLine = __LINE__;
                goto CBC_OPTIMIZATION_ERROR;
            }
            if ((model->isProvenInfeasible()) || (model->isProvenDualInfeasible())) {
                lp->status = LP_INTINFEASIBLE;
                goto CBC_OPTIMIZATION_CONCLUDED;
            }

            if (model->bestSolution())
                lp->status = LP_FEASIBLE;
            if (model->isProvenOptimal())
                lp->status = LP_OPTIMAL;

            if (model->bestSolution()) {
                memcpy(&((*(lp->_x))[0]), model->solver()->getColSolution(), sizeof(double)*lp_cols(lp));
                for ( int i=0 ; (i<lp_rows(lp)) ; ++i )
                {
                    double activity = model->getCbcRowActivity()[i];
                    double lower = model->getCbcRowLower()[i];
                    double upper = model->getCbcRowUpper()[i];
                    (*lp->_slack)[i] = std::min( upper-activity, activity-lower );
                }


                for (int i = 0; (i < model->numberSavedSolutions()) ; i++) {
                    const double *si = model->savedSolution(i);
                    vector< double > ti;
                    ti.clear();
                    ti.insert(ti.end(), si, si + lp_cols(lp));

                    lp->_savedSol->push_back(ti);
                    lp->_savedObj->push_back(model->savedSolutionObjective(i));
                } // saved solution
            } // best solution

            lp->obj = model->getObjValue();
        }
    }
CBC_OPTIMIZATION_CONCLUDED:
    if (deleteLP) {
        assert(linearProgram);
        delete linearProgram;
    }

    goto OPTIMIZATION_CONCLUDED;
CBC_OPTIMIZATION_ERROR:
    if (deleteLP) {
        assert(linearProgram);
        delete linearProgram;
    }

    goto OPTIMIZATION_ERROR;
#endif
#ifdef CPX
    bool restoreColumnTypes = false;
    vector< char > ctypes;

    //lp_mipstart_debug( lp );
    lp_config_cpx_params(lp);

    if ((isMIP) && (!lp->optAsContinuous)) {
        int cpxError = CPXmipopt(LPcpxDefaultEnv, lp->cpxLP);
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        int solStat = CPXgetstat(LPcpxDefaultEnv, lp->cpxLP);
        bool hasSolution = false;

        lp->status = LP_NO_SOL_FOUND;

        CPXsetdblparam( LPcpxDefaultEnv,  CPX_PARAM_EPAGAP, 1e-7 );
        CPXsetdblparam( LPcpxDefaultEnv,  CPX_PARAM_EPGAP, 1e-5 );

        switch (solStat) {
            case CPXMIP_OPTIMAL :
            case CPXMIP_OPTIMAL_TOL :
            {
                hasSolution = true;
                int status  = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
                if (status) {
                    fprintf(stderr, "Could not get objval. At %s:%d.\n", __FILE__, __LINE__);
                    abort();
                }

                lp->status = LP_OPTIMAL;
                break;
            }
            case CPX_STAT_INFEASIBLE :
            case CPX_STAT_INForUNBD :
            case CPXMIP_INFEASIBLE :
            case CPXMIP_INForUNBD :
                lp->status = LP_INFEASIBLE;
                break;
            case CPX_STAT_UNBOUNDED :
            case CPXMIP_UNBOUNDED :
                lp->status = LP_UNBOUNDED;
                break;
            default:
            {
                int status  = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
                if (status)
                    lp->status = LP_NO_SOL_FOUND;
                else {
                    hasSolution = true;
                    lp->status = LP_FEASIBLE;
                }
                break;
            }
        }

        if (hasSolution) {
            assert( ((int)lp->_x->size())>=lp_cols(lp) );
            // getting x
            int cpxError = CPXgetx(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_x))[0]), 0, lp_cols(lp) - 1);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
            // getting slack
            double *slack = &((*lp->_slack)[0]);
            cpxError = CPXgetslack(LPcpxDefaultEnv, lp->cpxLP, slack, 0, lp_rows(lp) - 1);
            for ( int i=0 ; (i<lp_rows(lp)) ; ++i )
                slack[i] = fabs(slack[i]);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

            /* getting solution pool */
            int nsols = CPXgetsolnpoolnumsolns( LPcpxDefaultEnv, lp->cpxLP );
            if (nsols)
            {
                lp->_savedSol->resize( nsols, vector< double >( lp_cols(lp), 0.0 )  );
                lp->_savedObj->resize( nsols, DBL_MAX );
                for ( int i=0 ; (i<nsols) ; ++i )
                {
                    int error = CPXgetsolnpoolobjval( LPcpxDefaultEnv, lp->cpxLP, i, &((*(lp->_savedObj))[i]) );
                    assert( !error );
                    cpxError = CPXgetsolnpoolx( LPcpxDefaultEnv, lp->cpxLP, i,  &((*(lp->_savedSol))[i][0]), 0, lp_cols(lp)-1 );
                    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
                }
            }
         }

        if ( ( lp->status != LP_INFEASIBLE && lp->status != LP_UNBOUNDED ) )
        {
            int cpxError = CPXgetbestobjval( LPcpxDefaultEnv, lp->cpxLP, &lp->bestBound );
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }

        goto OPTIMIZATION_CONCLUDED;
    }
    else {
        int status = 0;
        if (CPXgetprobtype(LPcpxDefaultEnv, lp->cpxLP) == CPXPROB_MILP) {
            /* backing up column types */
            ctypes.resize( lp_cols(lp), CPX_CONTINUOUS );
            int cpxError = CPXgetctype( LPcpxDefaultEnv, lp->cpxLP, &ctypes[0], 0, lp_cols(lp)-1 );
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
            /* changing problem type */
            cpxError = CPXchgprobtype(LPcpxDefaultEnv, lp->cpxLP, CPXPROB_LP);
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
            restoreColumnTypes = true;
        }
        int cpxError = CPXlpopt(LPcpxDefaultEnv, lp->cpxLP );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        int solStat = CPXgetstat(LPcpxDefaultEnv, lp->cpxLP);

        switch (solStat) {
            case CPX_STAT_OPTIMAL :
            {
                status = CPXgetobjval(LPcpxDefaultEnv, lp->cpxLP, &lp->obj);
                if (status) {
                    sprintf(errorMsg, "Could not get objval.");
                    errorLine = __LINE__;
                    goto OPTIMIZATION_ERROR;
                }

                int cpxError = CPXgetx(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_x))[0]), 0, lp_cols(lp) - 1);
                lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

                cpxError = CPXgetdj(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_rc))[0]), 0, lp_cols(lp) - 1);
                lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

                cpxError = CPXgetpi(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_pi))[0]), 0, lp_rows(lp) - 1);
                lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

                cpxError = CPXgetslack(LPcpxDefaultEnv, lp->cpxLP, &((*(lp->_slack))[0]), 0, lp_rows(lp) - 1);
                lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

                lp->status = LP_OPTIMAL;
                goto OPTIMIZATION_CONCLUDED;
            }

            case CPX_STAT_INFEASIBLE :
                                    lp->status = LP_INFEASIBLE;
                                    goto OPTIMIZATION_CONCLUDED;
                                    break;
            case CPX_STAT_INForUNBD :
                                    lp->status = LP_INFEASIBLE;
                                    goto OPTIMIZATION_CONCLUDED;
                                    break;
            case CPX_STAT_UNBOUNDED :
                                    lp->status = LP_UNBOUNDED;
                                    goto OPTIMIZATION_CONCLUDED;
                                    break;
            default :
                                    char statStr[256];
                                    CPXgetstatstring(LPcpxDefaultEnv, solStat, statStr);
                                    sprintf(errorMsg, "CPLEX CPXlpopt returned unhandled optimization status %s.\n", statStr);
                                    errorLine = __LINE__;
                                    goto OPTIMIZATION_ERROR;
                                    break;
        }
    }
#endif
#ifdef GRB
    lp->status = LP_ERROR;
    bool hasSolution = false;
    bool restoreColumnTypes = false;
    vector< char > origVTypes;

    int nIntVars = 0;
    int grbError = GRBgetintattr( lp->lp, "NumIntVars", &nIntVars );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    // setting up as a continuous model
    if ( lp->optAsContinuous && nIntVars )
    {
        // backing up
        origVTypes.resize( lp_cols(lp) );
        grbError = GRBgetcharattrarray( lp->lp, "VType", 0, lp_cols(lp), &origVTypes[0] );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        restoreColumnTypes = true;


        // changing all variables to continuous
        char *vTypes = new char[lp_cols(lp)];
        memset( vTypes, GRB_CONTINUOUS, sizeof(char)*lp_cols(lp) );
        grbError = GRBsetcharattrarray( lp->lp, "VType", 0, lp_cols(lp), vTypes );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        grbError = GRBupdatemodel( lp->lp );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        delete[] vTypes;
    }

    lp_config_grb_params( lp );
    grbError = GRBoptimize( lp->lp );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    int optStatus;
    grbError = GRBgetintattr(lp->lp, GRB_INT_ATTR_STATUS, &optStatus);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    switch (optStatus)
    {
        case GRB_OPTIMAL:
        {
            lp->status = LP_OPTIMAL;
            hasSolution = true;
            break;
        }
        case GRB_SUBOPTIMAL:
        {
            lp->status = LP_FEASIBLE;
            hasSolution = true;
            break;
        }
        case GRB_INFEASIBLE:
        {
            lp->status = LP_INFEASIBLE;
            goto OPTIMIZATION_CONCLUDED;
            break;
        }
        case GRB_INF_OR_UNBD:
        {
            lp->status = LP_INFEASIBLE;
            goto OPTIMIZATION_CONCLUDED;
            break;
        }
        case GRB_UNBOUNDED:
        {
            lp->status = LP_UNBOUNDED;
            goto OPTIMIZATION_CONCLUDED;
            break;
        }
        case GRB_NODE_LIMIT :
        {
            int nSol = 0, grbError = 0;
            grbError = GRBgetintattr( lp->lp, "SolCount", &nSol );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

            lp->status = (nSol>0) ? LP_FEASIBLE : LP_NO_SOL_FOUND;
            hasSolution = (nSol>0);
            break;
        }
        case GRB_TIME_LIMIT :
        {
            int nSol = 0, grbError = 0;
            grbError = GRBgetintattr( lp->lp, "SolCount", &nSol );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

            lp->status = (nSol>0) ? LP_FEASIBLE : LP_NO_SOL_FOUND;
            hasSolution = (nSol>0);
            break;
         }
        case GRB_SOLUTION_LIMIT :
        {
            int nSol = 0, grbError = 0;
            grbError = GRBgetintattr( lp->lp, "SolCount", &nSol );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
            lp->status = (nSol>0) ? LP_FEASIBLE : LP_NO_SOL_FOUND;
            hasSolution = (nSol>0);
            break;
        }
        case GRB_ITERATION_LIMIT:
        {
            int nSol = 0, grbError = 0;
            grbError = GRBgetintattr( lp->lp, "SolCount", &nSol );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
            lp->status = (nSol>0) ? LP_FEASIBLE : LP_NO_SOL_FOUND;
            hasSolution = (nSol>0);
            break;
        }
        case GRB_CUTOFF:
        {
            int nSol = 0, grbError = 0;
            grbError = GRBgetintattr( lp->lp, "SolCount", &nSol );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
            lp->status = (nSol>0) ? LP_FEASIBLE : LP_NO_SOL_FOUND;
            fprintf( stderr, "[warning]: specified cutoff is better than optimal solution. no solution information is available.\n");
            hasSolution = (nSol>0);
            break;
        }

         default:
            sprintf( errorMsg, "Unknow status: %d\n", optStatus );
            goto OPTIMIZATION_ERROR;
    }

    /* checking best bound for MIPs */
    if (nIntVars>0)
    {
        if (optStatus != GRB_OPTIMAL)
        {
            grbError = GRBgetdblattr( lp->lp, "ObjBound", &(lp->bestBound) );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        }
        else
        {
            grbError = GRBgetdblattr(lp->lp, "ObjVal", &lp->obj);
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
            lp->bestBound = lp->obj;
        }
    }

    if (hasSolution)
    {
        grbError = GRBgetdblattr(lp->lp, "ObjVal", &lp->obj);
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

        double *x = &((*lp->_x)[0]);
        grbError = GRBgetdblattrarray( lp->lp , GRB_DBL_ATTR_X, 0, lp_cols(lp), x );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

        double *slack = &((*lp->_slack)[0]);
        grbError = GRBgetdblattrarray( lp->lp , GRB_DBL_ATTR_SLACK, 0, lp_rows(lp), slack );
        for ( int i=0 ; (i<lp_rows(lp)) ; ++i )
            slack[i] = fabs(slack[i]);

        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

        // if it was solved as a linear program then dual information for rows should be available
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        if (nIntVars==0 || lp->optAsContinuous)
        {
            double *pi = &((*lp->_pi)[0]);
            grbError = GRBgetdblattrarray( lp->lp, GRB_DBL_ATTR_PI, 0, lp_rows(lp), pi );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

            double *rc = &((*lp->_rc)[0]);
            grbError = GRBgetdblattrarray( lp->lp , GRB_DBL_ATTR_RC, 0, lp_cols(lp), rc );
            lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        }

        goto OPTIMIZATION_CONCLUDED;
    }
#endif
OPTIMIZATION_CONCLUDED:

#ifdef CPX
    if (restoreColumnTypes)
    {
        vector< int > idx( lp_cols(lp), 0 );
        for ( int i=0 ; (i<lp_cols(lp)) ; i++ ) idx[i] = i;
        int cpxError = CPXchgctype( LPcpxDefaultEnv, lp->cpxLP, lp_cols(lp), &idx[0], &(ctypes[0]) );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
#endif
#ifdef GRB
    if (restoreColumnTypes)
    {
        grbError = GRBsetcharattrarray( lp->lp, "VType", 0, lp_cols(lp), &origVTypes[0] );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    }
#endif

    if ( (isMIP) && (lp->status == LP_OPTIMAL) )
        lp->bestBound = lp->obj;

    time_t endT; time(&endT);

    lp->solutionTime = difftime( endT, startT );
    if ((strlen(lp->solOutFN)) && (lp_obj_value(lp) != DBL_MAX))
        lp_write_sol(lp, lp->solOutFN);
    return lp->status;

OPTIMIZATION_ERROR:
    fprintf(stderr, "\n\n===--->>> ERROR <<<---===\n");
    fprintf(stderr, "\tAt lp.cpp, line: %d\n", errorLine);
    fprintf(stderr, "\tMessage: %s\n", errorMsg);
    fprintf(stderr, "\tSaving LP in error.lp\n");
    lp_write_lp(lp, "error.lp");
    abort();
}

double lp_obj_value(LinearProgram *lp)
{
    assert(lp != NULL);
    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    return lp->obj;
}

double *lp_row_price(LinearProgram *lp)
{
    assert(lp != NULL);

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    if (lp->status != LP_OPTIMAL) {
        fprintf(stderr, "\n\nERROR: no dual solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }

    return &((*(lp->_pi))[0]);
}

double *lp_row_slack( LinearProgram *lp )
{
    assert(lp != NULL);

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    if (lp->status != LP_OPTIMAL && lp->status!=LP_FEASIBLE) {
        fprintf(stderr, "\n\nERROR: no solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }

    return &((*(lp->_slack))[0]);
}

double *lp_x(LinearProgram *lp)
{
    assert(lp != NULL);

    if ((lp->status != LP_OPTIMAL) && (lp->status != LP_FEASIBLE)) {
        fprintf(stderr, "\n\nERROR: no solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    return &((*(lp->_x))[0]);
}

int lp_get_mip_emphasis(LinearProgram *lp)
{
    assert(lp != NULL);

    return lp->mipEmphasis;
}

void lp_set_mip_emphasis(LinearProgram *lp, const int mipEmphasis)
{
    assert(lp != NULL);

    lp->mipEmphasis = mipEmphasis;
}

void lp_printRootRelaxInfo(LinearProgram *lp)
{
    assert(lp != NULL);

    double obj = lp_obj_value(lp);
    printf("Root node linear relaxation info:\n");
    printf("\tObjective value: %g\n", obj);
    fflush(stdout);
}

#ifdef CBC
/* cut generation related includes */
#include <CglKnapsackCover.hpp>
#include <CglFlowCover.hpp>
#include <CglZeroHalf.hpp>
#include <CglMixedIntegerRounding.hpp>
#include <CglTwomir.hpp>
#include <CglLandP.hpp>
#include <CglRedSplit.hpp>
#include <CglGMI.hpp>

int lp_strengthen_with_cuts( LinearProgram *lp, const int maxRoundsCuts[] )
{
    OsiClpSolverInterface *osiLP = lp->_lp;
    int round = 1;
    int totalCuts = 0;
    assert( maxRoundsCuts );

    printf("optimizations %d obj %g\n", lp->nOptimizations, lp_obj_value(lp) );
#ifdef CBC
CUTGEN:
    int origRows = lp_rows( lp );
    int newCuts = 0;
    char message[256] = "";

    {
        if (!osiLP->optimalBasisIsAvailable())
            osiLP->initialSolve();

        OsiCuts cuts;

        if (round<=maxRoundsCuts[LPC_FLOW])
        {
            int oc = cuts.sizeCuts();
            CglFlowCover flowCoverCuts;
            flowCoverCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)FlowCover ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_ZERO_HALF])
        {
            int oc = cuts.sizeCuts();
            CglZeroHalf zhCuts;
            zhCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)ZeroHalf ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_MIR])
        {
            int oc = cuts.sizeCuts();
            CglMixedIntegerRounding mirCuts;
            mirCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)MIR ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_GOMORY])
        {
            int oc = cuts.sizeCuts();
            CglGMI gomoryCuts;
            gomoryCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)Gomory ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_KNAPSACK])
        {
            int oc = cuts.sizeCuts();
            CglKnapsackCover coverCuts;
            coverCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)Knapsack ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_TWO_MIR])
        {
            int oc = cuts.sizeCuts();
            CglTwomir twoMirCuts;
            twoMirCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)TwoMIR ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_L_AND_P])
        {
            int oc = cuts.sizeCuts();
            CglLandP landPCuts;
            landPCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)LiftAndProject ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        if (round<=maxRoundsCuts[LPC_REDUCE])
        {
            int oc = cuts.sizeCuts();
            CglRedSplit reduceCuts;
            reduceCuts.generateCuts( *osiLP, cuts );
            if (cuts.sizeCuts()>oc)
            {
                char str[64];
                sprintf( str, "(%d)ReduceAndSplit ", cuts.sizeCuts()-oc );
                strcat( message, str );
            }
        }

        newCuts = cuts.sizeCuts();
        osiLP->applyCuts( cuts );
        assert( newCuts == lp_rows(lp)-origRows );
    }


    if (newCuts)
    {
        if (!lp->silent)
            printf( "%d new cuts [%s] inserted in formulation.\n", newCuts, message );
        totalCuts += newCuts;
        osiLP->resolve();
        if (!lp->silent)
            printf( "objective value is now %g\n\n", lp_obj_value(lp) );
        ++round;
        goto CUTGEN;
    }
    else
    {
        if (!lp->silent)
            printf( "no violated cuts found.\n" );
    }

    return totalCuts;
#endif
    printf("WARNING:  this solver does not supports the generation of cuts outside the branch and cut tree.\n");
}

void lp_set_callback( LinearProgram *lp, lp_cb callback, void *data )
{
    lp->callback_ = callback;
    lp->data_ = data;
}
#endif // CBC

void lp_free(LinearProgramPtr *lp)
{
    assert(lp != NULL);
    assert(*lp != NULL);

#ifdef GLPK
    glp_delete_prob((*lp)->_lp);
#endif
#ifdef CBC
    if ((*lp)->justWrapOsi == 0) {
        if ((*lp)->cglPP == NULL) {
            delete (*lp)->_lp;
            (*lp)->osiLP = NULL;
        }
        else {
            delete (*lp)->cglPP;
            (*lp)->_lp = NULL;
            (*lp)->osiLP = NULL;
        }
    }

    if ( ((*lp)->cutPool)&&((*lp)->ownCutPool) )
        delete ((*lp)->cutPool);
#endif
#ifdef CPX
    //CPXfree((*lp)->cpxLP);//LPcpxDefaultEnv);
    //CPXfreeparenv(LPcpxDefaultEnv,  &LPcpxDefaultEnv );
    CPXfreeprob(LPcpxDefaultEnv, &((*lp)->cpxLP));
#endif // CPX
#ifdef GRB
  int grbError = GRBfreemodel( ((*lp)->lp) );
  lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif

    delete (*lp)->_x;
    delete (*lp)->_rc;
    delete (*lp)->_obj;
    delete (*lp)->_idx;
    delete (*lp)->_coef;
    delete (*lp)->_pi;
    delete (*lp)->_slack;
    delete (*lp)->_mipStart;
    delete (*lp)->_savedSol;
    delete (*lp)->_savedObj;
    delete (*lp)->colNameIdx;
    delete (*lp)->_orig;
    delete (*lp)->msIdx;
    delete (*lp)->msVal;

    free(*lp);
    *lp = NULL;
}

int lp_col( LinearProgram *lp, int col, int *idx, double *coef )
{
    LP_CHECK_COL_INDEX( lp, col );

    int result = -INT_MAX;

#ifdef CPX
    int surplus = -INT_MAX;
    int rmatbeg[] = { 0, lp_rows(lp) };
    int cpxError = CPXgetcols(LPcpxDefaultEnv, lp->cpxLP, &result, &rmatbeg[0], idx, coef, lp_rows(lp), &surplus, col, col );
    assert( surplus>=0 );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef GLPK
    result = glp_get_mat_col(lp->_lp, col + 1, idx - 1, coef - 1);
    for ( int i = 0 ; (i<result) ; ++i )
        idx[i]--;
#endif
#ifdef CBC
    const CoinPackedMatrix *cpmCol =  lp->osiLP->getMatrixByCol();
    const int nzCol = cpmCol->getVectorLengths()[col];
    const CoinBigIndex *starts = cpmCol->getVectorStarts();
    const int *ridx = cpmCol->getIndices() + starts[col];
    const double *rcoef = cpmCol->getElements() + starts[col];

    for (int j = 0 ; (j < nzCol) ; ++j) {
        idx[j] = ridx[j];
        coef[j] = rcoef[j];
    }

    result = nzCol;
#endif
#ifdef GRB
    int nz = 0;
    int beg[2];
    int grbError = GRBgetvars( lp->lp, &nz, beg,  idx, coef, col, 1 );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    result = nz;
#endif

    return result;
}

int lp_row(LinearProgram *lp, int row, int *idx, double *coef)
{
    LP_CHECK_ROW_INDEX( lp, row );

    int result = -INT_MAX;

#ifdef GRB
    int nz = 0;
    int beg[2];
    int grbError = GRBgetconstrs( lp->lp, &nz, beg,  idx, coef, row, 1 );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    return nz;
#endif
#ifdef CPX
    int surplus= -INT_MAX;
    int rmatbeg[2] = { -INT_MAX, -INT_MAX };
    int cpxError = CPXgetrows(LPcpxDefaultEnv, lp->cpxLP, &result, &rmatbeg[0], idx, coef, lp_cols(lp) + 1, &surplus, row, row);
    assert( surplus>=0 );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif

#ifdef CBC
    const CoinPackedMatrix *cpmRow =  lp->osiLP->getMatrixByRow();
    const int nzRow = cpmRow->getVectorLengths()[row];
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const int *ridx = cpmRow->getIndices() + starts[row];
    const double *rcoef = cpmRow->getElements() + starts[row];
    for (int j = 0 ; (j < nzRow) ; ++j) {
        idx[j] = ridx[j];
        coef[j] = rcoef[j];
    }

    result = nzRow;
#endif
#ifdef GLPK
    result = glp_get_mat_row(lp->_lp, row + 1, idx - 1, coef - 1);
    for (int i = 0 ; (i<result) ; ++i)
        idx[i]--;
#endif

    return result;
}

double lp_rhs(LinearProgram *lp, int row)
{
    assert(lp != NULL);
    LP_CHECK_ROW_INDEX( lp, row );

#ifdef GRB
    double rhs;
    int grbError = GRBgetdblattrelement( lp->lp, "RHS", row,  &rhs );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    return rhs;
#endif
#ifdef DEBUG
    assert(row >= 0);
    assert(row < lp_rows(lp));
#endif
#ifdef CPX
    double rhs;
    int cpxError = CPXgetrhs(LPcpxDefaultEnv, lp->cpxLP, &rhs, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return rhs;
#endif
#ifdef CBC
    return lp->osiLP->getRightHandSide()[row];
#endif
#ifdef GLPK
    switch (glp_get_row_type(lp->_lp, row + 1)) {
        case GLP_LO:
            return glp_get_row_lb(lp->_lp, row + 1);
        case GLP_UP:
            return glp_get_row_ub(lp->_lp, row + 1);
        case GLP_FX:
            return glp_get_row_ub(lp->_lp, row + 1);
        default:
            abort();
    }
#endif
}

char lp_sense(LinearProgram *lp, int row)
{
    assert(lp != NULL);
    LP_CHECK_ROW_INDEX( lp, row );

#ifdef GRB
    char sense;
    int grbError = GRBgetcharattrelement( lp->lp, "Sense", row, &sense );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    switch (sense)
    {
        case '<':
            return 'L';
        case '>':
            return 'G';
        case '=':
            return 'E';
        default:
            fprintf( stderr, "Sense %c not recognized.\n", sense );
            exit( EXIT_FAILURE );
    }
#endif
#ifdef DEBUG
    assert(row >= 0);
    assert(row < lp_rows(lp));
#endif
#ifdef CPX
    int cpxError;
    char result;
    cpxError = CPXgetsense(LPcpxDefaultEnv, lp->cpxLP, &result, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return result;
#endif
#ifdef CBC
    return lp->osiLP->getRowSense()[row];
#endif
#ifdef GLPK
    switch (glp_get_row_type(lp->_lp, row + 1)) {
        case GLP_LO:
            return 'G';
        case GLP_UP:
            return 'L';
        case GLP_FX:
            return 'E';
        default:
            abort();
    }
#endif
}

char *lp_row_name(LinearProgram *lp, int row, char *dest)
{
    assert(lp != NULL);

    LP_CHECK_ROW_INDEX(lp, row);
    assert(dest);

#ifdef GRB
    char *name;
    int grbError =
        GRBgetstrattrelement( lp->lp, "ConstrName", row, &name);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    strcpy( dest, name );
    return dest;
#endif
#ifdef CPX
    int surplus = 0;
    int cpxError = CPXgetrowname(LPcpxDefaultEnv, lp->cpxLP, &dest, dest, 256, &surplus, row, row);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    strcpy(dest, lp->osiLP->getRowName(row).c_str());
#endif
#ifdef GLPK
    strcpy(dest, glp_get_row_name(lp->_lp, row + 1));
#endif
    return dest;
}

char *lp_col_name(LinearProgram *lp, int col, char *dest)
{
    assert(lp != NULL);

    LP_CHECK_COL_INDEX(lp, col);
    assert(dest);

#ifdef GRB
    char *name;
    int grbError =
        GRBgetstrattrelement( lp->lp, GRB_STR_ATTR_VARNAME, col, &name);
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    strcpy( dest, name );

    return dest;
#endif
#ifdef CPX
    int surplus = 0;
    int cpxError = CPXgetcolname(LPcpxDefaultEnv, lp->cpxLP, &dest, dest, 256, &surplus, col, col);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    strcpy(dest, lp->osiLP->getColName(col).c_str());
#endif
#ifdef GLPK
    strcpy(dest, glp_get_col_name(lp->_lp, col + 1));
#endif

    return dest;
}

double lp_col_lb(LinearProgram *lp, int col)
{
    assert(lp != NULL);
    LP_CHECK_COL_INDEX(lp, col);

#ifdef GRB
    double lb;
    int grbError = GRBgetdblattrelement( lp->lp, GRB_DBL_ATTR_LB, col, &lb );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return lb;
#endif
#ifdef DEBUG
    assert(col >= 0);
    assert(col < lp_cols(lp));
#endif
#ifdef CPX
    double lb;
    int begin = col;
    int end = col;
    int cpxError = CPXgetlb(LPcpxDefaultEnv, lp->cpxLP, &lb, begin, end);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    return lb;
#endif
#ifdef CBC
    return lp->osiLP->getColLower()[col];
#endif
#ifdef GLPK
    return glp_get_col_lb(lp->_lp, col + 1);
#endif

}

double lp_col_ub(LinearProgram *lp, int col)
{
    assert(lp != NULL);
    LP_CHECK_COL_INDEX(lp, col);

#ifdef GRB
    double ub;
    int grbError = GRBgetdblattrelement( lp->lp, GRB_DBL_ATTR_UB, col, &ub );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return ub;
#endif
#ifdef DEBUG
    assert(col >= 0);
    assert(col < lp_cols(lp));
#endif
#ifdef CPX
    double ub;
    int begin = col;
    int end = col;
    int cpxError = CPXgetub(LPcpxDefaultEnv, lp->cpxLP, &ub, begin, end);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    return ub;
#endif
#ifdef CBC
    return lp->osiLP->getColUpper()[col];
#endif
#ifdef GLPK
    return glp_get_col_ub(lp->_lp, col + 1);
#endif
}

void lp_set_obj(LinearProgram *lp, double *obj)
{
    assert(lp != NULL);

#ifdef GRB
    int grbError = GRBsetdblattrarray( lp->lp, "Obj", 0, lp_cols(lp) , obj );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    grbError = GRBupdatemodel( lp->lp );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif
#ifdef CPX
    vector< int > idx(lp_cols(lp), 0);
    for (int i = 0 ; (i < lp_cols(lp)) ; i++)
        idx[i] = i;

    int cpxError = CPXchgobj(LPcpxDefaultEnv, lp->cpxLP, lp_cols(lp), &idx[0], obj);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    return lp->osiLP->setObjective(obj);
#endif
#ifdef GLPK
    for (int i = 0 ; (i < lp_cols(lp)) ; ++i)
        glp_set_obj_coef(lp->_lp, i + 1, obj[i]);
#endif
}

void lp_add_col(LinearProgram *lp, double obj, double lb, double ub, char integer, char *name, int nz, int *rowIdx, double *rowCoef)
{
    assert(lp != NULL);

#ifdef GRB
    char vType = GRB_CONTINUOUS;
    if (integer)
    {
        if ( abs( lb ) < 1e-10 && abs(ub-1.0) < 1e-10 )
            vType = GRB_BINARY;
        else
            vType = GRB_INTEGER;
    }
    if (ub == DBL_MAX)
        ub = GRB_INFINITY;
    if ( lb == -DBL_MAX )
        lb = -GRB_INFINITY;

    int grbError = GRBaddvar( lp->lp, nz, rowIdx, rowCoef, obj, lb, ub, vType, name );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    lp->tmpCols++;
#endif
#ifdef CPX
    char type = CPX_CONTINUOUS;
    if (integer) {
        if ((fabs(lb) < EPS) && (fabs(ub - 1.0) < EPS))
            type = CPX_BINARY;
        else
            type = CPX_INTEGER;
    }

    int matBeg[] = { 0, nz };
    int cpxError = CPXaddcols( LPcpxDefaultEnv, lp->cpxLP, 1, nz, &obj, matBeg, rowIdx, rowCoef, &lb, &ub, &name );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);

    if ( type != CPX_CONTINUOUS )
    {
        int idx = lp_cols(lp)-1;
        cpxError = CPXchgctype( LPcpxDefaultEnv, lp->cpxLP, 1, &idx, &type  );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }

    return;
#endif // CPX
#ifdef CBC
    int starts[] = { 0, nz };
    lp->osiLP->addCols(1, starts, rowIdx, rowCoef, &lb, &ub, &obj);
    if (integer)
        lp->osiLP->setInteger(lp_cols(lp) - 1);
    lp->osiLP->setColName(lp_cols(lp) - 1, name);
#endif
#ifdef GLPK
    int cols;

    glp_add_cols(lp->_lp, 1);
    cols = lp_cols(lp);

    if (name)
        glp_set_col_name(lp->_lp, cols, name);

    glp_set_obj_coef(lp->_lp, cols, obj);

    if (integer) {
        if ((fabs(ub - 1.0) <= EPS) && (fabs(lb) <= EPS))
            glp_set_col_kind(lp->_lp, cols, GLP_BV);
        else
            glp_set_col_kind(lp->_lp, cols, GLP_IV);
    }

    if ((lb != -DBL_MAX) && (ub != DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_DB, lb , ub);
    else if ((ub == DBL_MAX) && (lb != -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_LO, lb , ub);
    else if ((ub != DBL_MAX) && (lb == -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_UP, lb , ub);
    else if ((ub == DBL_MAX) && (lb == -DBL_MAX))
        glp_set_col_bnds(lp->_lp, cols, GLP_FR, -DBL_MAX , DBL_MAX);

    for (int i = 0 ; (i < nz) ; ++i)
        rowIdx[i]++;

    glp_set_mat_col(lp->_lp, cols, nz, rowIdx - 1, rowCoef - 1);

    for (int i = 0 ; (i < nz) ; ++i)
        rowIdx[i]--;
#endif

#ifdef GRB
    int idxCol = lp_cols(lp) +lp->tmpCols;
#else
    int idxCol = lp_cols(lp);
#endif
    map< string, int > &mNames = (*lp->colNameIdx);
    if ( mNames.find( name ) != mNames.end() )
    {
        fprintf( stderr, "ERROR: adding variable with repeated name %s\n", name );
        abort();
    }
    mNames[name] = idxCol;
}

void lp_set_rhs(LinearProgram *lp, int row, double rhs)
{
    assert(lp != NULL);
    LP_CHECK_ROW_INDEX( lp, row );

#ifdef GRB
    int grbError = GRBsetdblattrelement( lp->lp, GRB_DBL_ATTR_RHS, row, rhs );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    lp->nModelChanges++;
    return;
#endif

#ifdef CPX
    int cpxError = CPXchgrhs(LPcpxDefaultEnv, lp->cpxLP, 1, &row, &rhs);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#else


#ifdef CBC
    char sense = lp_sense(lp, row);
    switch (sense) {
        case 'E':
            lp->osiLP->setRowBounds(row, rhs, rhs);
            break;
        case 'G':
            lp->osiLP->setRowBounds(row, rhs, COIN_DBL_MAX);
            break;
        case 'L':
            lp->osiLP->setRowBounds(row, -COIN_DBL_MAX, rhs);
            break;
        default:
            fprintf(stderr, "Unknow sense: %c!\n", sense);
            abort();
            exit(1);
    }
#endif
#ifdef GLPK
    char sense = lp_sense(lp, row);
    switch (sense) {
        case 'E':
            glp_set_row_bnds(lp->_lp, row + 1, GLP_FX, rhs, rhs);
            break;
        case 'G':
            glp_set_row_bnds(lp->_lp, row + 1, GLP_LO, rhs, DBL_MAX);
            break;
        case 'L':
            glp_set_row_bnds(lp->_lp, row + 1, GLP_UP, -DBL_MAX, rhs);
            break;
        default :
            fprintf(stderr, "Unknow sense: %c!\n", sense);
            abort();
            exit(1);
    }
#endif
#endif
}

double lp_solution_time(LinearProgram *lp)
{
    assert(lp != NULL);

    return lp->solutionTime;
}

#ifdef CPX
void lp_check_for_cpx_error(CPXENVptr env, int errorCode, const char *sourceFile, int sourceLine)
{
    if (errorCode) {
        char errorStr[256];
        CPXgeterrorstring(LPcpxDefaultEnv, errorCode, errorStr);
        fprintf(stderr, "CPLEX Error: %s\n", errorStr);
        fprintf(stderr, "Inside LP Library - %s:%d\n\n", sourceFile, sourceLine);
        abort();
    }
}
#endif

#ifdef GRB
void lp_check_for_grb_error(GRBenv* env, int errorCode, const char *sourceFile, int sourceLine)
{
    if (errorCode) {
        fprintf(stderr, "Gurobi Error: %s\n", GRBgeterrormsg( LPgrbDefaultEnv ) );
        fprintf(stderr, "Inside LP Library - %s:%d\n\n", sourceFile, sourceLine);
        abort();
    }
}
#endif

void lp_set_cuts(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->cuts = onOff;
}

void lp_set_max_seconds(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSeconds = _max;
}

void lp_set_max_solutions(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSolutions = _max;
}


int lp_num_saved_sols( LinearProgram *lp )
{
    assert(lp != NULL);
    return (int)lp->_savedSol->size();

}

void lp_set_max_saved_sols(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxSavedSols = _max;
}

void lp_set_max_nodes(LinearProgram *lp, int _max)
{
    assert(lp != NULL);
    lp->maxNodes = _max;
}

void lp_set_print_messages(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->printMessages = onOff;
}

void lp_set_parallel(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->parallel = onOff;
}

void lp_set_heur_proximity(LinearProgram *lp, char onOff)
{
    assert(lp != NULL);
    lp->heurProx = onOff;
}

void lp_set_heur_fp_passes(LinearProgram *lp, int passes)
{
    assert(lp != NULL);

    lp->heurFPPasses = passes;
}

#ifdef CBC
void lp_config_cbc_params(LinearProgram *lp, vector<string> &cbcP)
{
    assert(lp != NULL);

    cbcP.push_back("someMIP"); // problem name
    if (lp->cuts != INT_NOT_SET) {
        if (lp->cuts)
            cbcP.push_back("-cuts=on");
        else
            cbcP.push_back("-cuts=off");
    }
    if (lp->printMessages != INT_NOT_SET) {
        if (lp->printMessages)
            cbcP.push_back("-log=1");
        else
            cbcP.push_back("-log=0");
    }

    cbcP.push_back("-zero=ifmove");
    cbcP.push_back("-multiple=2");

    if ( lp->parallel != INT_NOT_SET )
    {
        if ( lp->parallel )
        {
#ifdef _OPENMP
            int nthreads = omp_get_num_procs();
#else
            int nthreads = 4;
#endif
#ifdef CBC
            if (nthreads>4)
                nthreads = 4; // to avoid excessive memory use
#endif
            if (lp->printMessages!=0)
                printf("CBC will use %d threads.\n", nthreads );
            char stropt[256];  sprintf( stropt, "-threads=%d", nthreads );
            cbcP.push_back( stropt );
        }
        else
            cbcP.push_back( "-threads=0" );
    }

    if (lp->maxSeconds != INT_NOT_SET)
    {
        cbcP.push_back("-timeM=elapsed");
        cbcP.push_back("-seconds=" + SSTR(lp->maxSeconds));
    }
    if (lp->maxSolutions != INT_NOT_SET)
        cbcP.push_back("-maxSol=" + SSTR(lp->maxSolutions));
    if (lp->maxNodes != INT_NOT_SET)
        cbcP.push_back("-maxNodes=" + SSTR(lp->maxNodes));
    if (lp->heurFPPasses  != INT_NOT_SET)
        cbcP.push_back("-passF=" + SSTR(lp->heurFPPasses));
    if (lp->heurProx != INT_NOT_SET) {
        if (lp->heurProx)
            cbcP.push_back("-proximity=on");
        else
            cbcP.push_back("-proximity=off");
    }
    if (lp->maxSavedSols != INT_NOT_SET)
        cbcP.push_back("-maxSaved=" + SSTR(lp->maxSavedSols));

    if (strlen(lp->solInFN))
        cbcP.push_back("-mips=" + string(lp->solInFN));

    if (lp->absMIPGap!=DBL_MAX)
    {
        char pstr[256];
        sprintf( pstr, "-allowableGap=%g", lp->absMIPGap );
        cbcP.push_back(pstr);
    }
    if (lp->relMIPGap!=DBL_MAX)
    {
        char pstr[256];
        sprintf( pstr, "-ratioGap=%g", lp->relMIPGap );
        cbcP.push_back(pstr);
    }
    if ( lp->mipPreProcess != CHAR_MAX )
    {
        if (!lp->mipPreProcess)
            cbcP.push_back("-preprocess=off");
    }

    cbcP.push_back("-solve");
    cbcP.push_back("-quit");
}
#endif

void lp_set_col_bounds(LinearProgram *lp, int col, const double lb, const double ub)
{
    LP_CHECK_COL_INDEX(lp, col);

    double l = lb;
    double u = ub;
    if ( lp_is_integer(lp,col) )
    {
        if (( l!=-DBL_MAX ) && (l!=DBL_MIN))
            l = floor(lb + 0.5);
        if ( u!=DBL_MAX )
            u = floor(ub + 0.5);
    }

#ifdef GRB
    int grbStatus = 0;
    grbStatus = GRBsetdblattrelement( lp->lp, GRB_DBL_ATTR_LB, col, lb );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbStatus, __FILE__, __LINE__ );
    grbStatus = GRBsetdblattrelement( lp->lp, GRB_DBL_ATTR_UB, col, ub );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbStatus, __FILE__, __LINE__ );
    lp->nModelChanges++;
#endif
#ifdef CBC
    OsiSolverInterface *linearProgram = lp->osiLP;
    linearProgram->setColBounds(col, l, u);
#endif
#ifdef GLPK
    if ( fabs(l-u)<= EPS )
        glp_set_col_bnds(lp->_lp, col+1, GLP_FX, l, l);
    else
        if ( (l==-DBL_MAX) || (l==DBL_MIN) )
            glp_set_col_bnds( lp->_lp, col+1, GLP_UP, l, u );
        else
            if ( u==DBL_MAX )
                glp_set_col_bnds( lp->_lp, col+1, GLP_LO, l, u );
            else
                glp_set_col_bnds( lp->_lp, col+1, GLP_DB, l, u );
#endif
#ifdef CPX
    const int idx[] = { col };
    if ( fabs(l-u)<= EPS )
    {
        int cpxError = CPXchgbds( LPcpxDefaultEnv, lp->cpxLP, 1, idx, "B", &lb ) ;
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
    else
    {
        if ( (l!=-DBL_MAX) && (l!=DBL_MIN) )
        {
            int cpxError;
            cpxError = CPXchgbds( LPcpxDefaultEnv, lp->cpxLP, 1, idx, "L", &lb ) ;
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }
        if ( u!=DBL_MAX )
        {
            int cpxError = CPXchgbds( LPcpxDefaultEnv, lp->cpxLP, 1, idx, "U", &ub ) ;
            lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
        }
    }
#endif


}

double lp_saved_sol_obj(LinearProgram *lp, int isol)
{
#ifdef GRB
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif
    assert(isol >= 0);
    assert(isol < (int)lp->_savedObj->size());
    return lp->_savedObj->at(isol);
}

double *lp_saved_sol_x(LinearProgram *lp, int isol)
{
#ifdef GRB
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif
    assert(lp != NULL);

    assert(isol >= 0);
    assert(isol < (int)lp->_savedSol->size());
    return &(lp->_savedSol->at(isol)[0]);
}

LinearProgram *lp_clone(LinearProgram *lp)
{
    assert(lp != NULL);

    LinearProgram *result = lp_create();

    memcpy(result, lp, sizeof(LinearProgram));

    lp_initialize(result);

    (*result->_x)         = (*lp->_x);
    (*result->_rc)        = (*lp->_rc);
    (*result->_obj)       = (*lp->_obj);
    (*result->_pi)        = (*lp->_pi);
    (*result->_slack)     = (*lp->_slack);
    (*result->_mipStart)  = (*lp->_mipStart);
    (*result->_savedSol)  = (*lp->_savedSol);
    (*result->_savedObj)  = (*lp->_savedObj);
    (*result->colNameIdx) = (*lp->colNameIdx);
    (*result->_orig)      = (*lp->_orig);
    result->callback_     = lp->callback_;
    result->data_         = lp->data_;

#ifdef GRB
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif

#ifdef CBC
    result->osiLP = lp->osiLP->clone();
    result->_lp = dynamic_cast<OsiClpSolverInterface *>(result->osiLP);
#endif
#ifdef GLPK
    result->_lp = glp_create_prob();
    glp_copy_prob(result->_lp, lp->_lp, GLP_ON);
#endif
#ifdef CPX
    int cpxError;
    result->cpxLP = CPXcreateprob(LPcpxDefaultEnv, &cpxError, "mip");
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    //cpxError = CPXcopylp( LPcpxDefaultEnv, lp->cpxLP, lp_cols(lp), lp_rows(lp), CPXgetobjsen( lp->cpxLP), , const double * rhs, const char * sense, const int * matbeg, const int * matcnt, const int * matind, const double * matval, const double * lb, const double * ub, const double * rngval)
#endif

    return result;
}

void lp_fix_col(LinearProgram *lp, int col, double val)
{
    LP_CHECK_COL_INDEX(lp, col);

    if (lp_is_integer(lp, col))
        val = floor(val + 0.5);
#ifdef GRB
    int grbError = GRBsetdblattrelement( lp->lp, "LB", col, val );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    grbError = GRBsetdblattrelement( lp->lp, "UB", col, val );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    lp->nModelChanges++;
#endif
#ifdef CBC
    lp->osiLP->setColBounds(col, val, val);
#endif
#ifdef CPX
    char lu = 'B';
    int cpxError = CPXchgbds(LPcpxDefaultEnv , lp->cpxLP, 1, &col, &lu, &val);
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef GLPK
    glp_set_col_bnds(lp->_lp, col+1, GLP_FX, val, val);
#endif
}

LinearProgram *lp_pre_process(LinearProgram *lp)
{
    assert(lp != NULL);

#ifdef GRB
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif
#ifdef CBC
    LinearProgram *result = (LinearProgramPtr) malloc(sizeof(LinearProgram));

    lp_initialize(result);

    result->cglPP = new CglPreProcess();
    result->_lp = dynamic_cast<OsiClpSolverInterface *>(result->cglPP->preProcess(*(lp->_lp), false, 4));
    result->osiLP = dynamic_cast<OsiSolverInterface *>(result->_lp);
    result->_lp->setIntParam(OsiNameDiscipline, 1);
    result->_lp->messageHandler()->setLogLevel(0);
    result->_orig->resize( lp_cols(result), INT_MAX );
    memcpy( &((*(result->_orig))[0]), (result->cglPP->originalColumns()), sizeof(int)*lp_cols(result) );

    return result;
#else
    fprintf(stderr, "\nERROR: Pre-processing not implemented for other solvers.\n");
    abort();
#endif
}

void lp_initialize(LinearProgram *lp)
{
    assert(lp != NULL);

    lp->_x = new vector< double >();
    lp->_rc = new vector< double >();
    lp->_obj = new vector< double >();
    lp->_coef = new vector< double >();
    lp->_idx = new vector< int >();
    lp->_pi = new vector< double >();
    lp->_slack =  new vector< double >();
    lp->_mipStart = new std::vector< std::pair< std::string, double > >();
    lp->_savedSol = new vector< vector<double> >();
    lp->_savedObj = new vector<double>();
    lp->colNameIdx = new map< string, int >();
    lp->_orig = new vector< int >();
    lp->msIdx = new vector< int >();
    lp->msVal = new vector< double >();
    lp->callback_ = NULL;
    lp->data_ = NULL;

#ifdef CBC
    lp->cutPool = NULL;
    lp->ownCutPool = 0;
    lp->justWrapOsi = 0;
#endif

    lp->optAsContinuous = 0;
    lp->mipPreProcess = CHAR_MAX;
    lp->nOptimizations = 0;
    lp->silent = 0;
    lp->solutionTime = 0.0;
    lp->obj = DBL_MAX;
    lp->status = LP_ERROR;
    lp->absMIPGap = DBL_MAX;
    lp->relMIPGap = DBL_MAX;

    strcpy(lp->solOutFN, "");
    strcpy(lp->solInFN, "");

    /* parameters */
    lp->maxSeconds    = INT_NOT_SET;
    lp->maxSavedSols  = INT_NOT_SET;
    lp->heurFPPasses  = INT_NOT_SET;
    lp->heurProx      = INT_NOT_SET;
    lp->maxNodes      = INT_NOT_SET;
    lp->cuts          = INT_NOT_SET;
    lp->printMessages = INT_NOT_SET;
    lp->maxSolutions  = INT_NOT_SET;
    lp->mipEmphasis   = LP_ME_DEFAULT;
    lp->parallel      = INT_NOT_SET;
}

int lp_col_index(LinearProgram *lp, const char *name)
{
    assert(lp != NULL);

    map< string, int >::const_iterator mIt;
    mIt = lp->colNameIdx->find(string(name));
    if (mIt == lp->colNameIdx->end())
        return -1;

    return mIt->second;
}

#ifdef GLPK
void lp_config_glpk_params(LinearProgram *lp, glp_iocp *iocp)
{
    glp_init_iocp(iocp);
    if (lp->silent)
        iocp->msg_lev = GLP_MSG_OFF;

    if ( lp->cuts == INT_NOT_SET )
        lp->cuts = 1;

    if ( lp->cuts )
    {
        iocp->gmi_cuts  = GLP_ON;
        iocp->mir_cuts  = GLP_ON;
        iocp->cov_cuts  = GLP_ON;
        iocp->clq_cuts  = GLP_ON;
    }

    if ( lp->heurProx != INT_NOT_SET )
        if ( lp->heurProx )
            iocp->ps_heur = GLP_ON;
    if ( lp->heurFPPasses != INT_NOT_SET )
        if ( lp->heurFPPasses>=1 )
            iocp->fp_heur = GLP_ON;

    /*iocp->ps_heur = GLP_ON;*/
    iocp->presolve = GLP_ON;
    iocp->br_tech = GLP_BR_PCH;

    if ( lp->mipEmphasis != LP_ME_DEFAULT )
    {
        switch (lp->mipEmphasis) {
            case LP_ME_OPTIMALITY:
                iocp->presolve  = GLP_ON;
                iocp->pp_tech   = GLP_PP_ALL;
                iocp->gmi_cuts  = GLP_ON;
                iocp->mir_cuts  = GLP_ON;
                iocp->cov_cuts  = GLP_ON;
                iocp->clq_cuts  = GLP_ON;
                iocp->fp_heur   = GLP_ON;

                break;
        }
    }
    if (lp->maxSeconds != INT_NOT_SET)
        iocp->tm_lim = lp->maxSeconds * 1000.0;

}
#endif


#ifdef GRB
void lp_config_grb_params(LinearProgram *lp)
{
    assert( lp && lp->lp );

    GRBenv *env = GRBgetenv( lp->lp );
    if ( lp->maxSeconds != INT_NOT_SET )
        GRBsetdblparam( env, GRB_DBL_PAR_TIMELIMIT, lp->maxSeconds);

    if ( lp->maxSolutions != INT_NOT_SET )
        GRBsetintparam( env, GRB_INT_PAR_SOLUTIONLIMIT, lp->maxSolutions );

    if ( lp->maxNodes != INT_NOT_SET )
        GRBsetdblparam( env, GRB_DBL_PAR_NODELIMIT, lp->maxNodes );

    if ((lp->silent) || (!lp->printMessages))
        GRBsetintparam( env, GRB_INT_PAR_OUTPUTFLAG, 0 );
    else
        GRBsetintparam( env, GRB_INT_PAR_OUTPUTFLAG, 1 );

    if ( lp->absMIPGap != DBL_MAX )
        GRBsetdblparam( env, GRB_DBL_PAR_MIPGAPABS, lp->absMIPGap );

    if ( lp->relMIPGap != DBL_MAX )
        GRBsetdblparam( env, GRB_DBL_PAR_MIPGAP, lp->relMIPGap );

    if ( strlen( lp->solInFN ) )
    {
#define STR_SIZE 256
        vector< pair< string, double > > &mipStart = *(lp->_mipStart);
        char **colNames = (char**) malloc( sizeof(char*)*lp_cols(lp) );
        colNames[0] = (char*) calloc( lp_cols(lp)*STR_SIZE, sizeof(char) );
        for ( int i=1 ; i<lp_cols(lp) ; ++i )
            colNames[i] = colNames[i-1]+STR_SIZE;
        double *colValues = (double*) calloc( lp_cols(lp), sizeof(double) );

        int rCols = lp_read_mip_start( lp->solInFN, colNames, colValues );
        mipStart.clear();
        for ( int i=0 ; i<rCols ; ++i )
            mipStart.push_back( pair< string, double >( colNames[i], colValues[i] ) );

        free( colValues );
        free( colNames[0] );
        free( colNames );
#undef STR_SIZE
        lp->msIdx->clear(); lp->msVal->clear();
    }
    if ( (lp->_mipStart) && (lp->_mipStart->size()) )
        lp_gen_mipstart(lp);

    if ( (lp->msIdx->size()) )
    {
        printf("setting grb mips start\n");

        int *idx = &((*(lp->msIdx))[0]);
        double *coef = &((*(lp->msVal))[0]);

        const int nNz = lp->msIdx->size();
        /* checking contents of msIdx */
        for ( int i=0 ; i<nNz ; ++i )
        {
            assert( idx[i]>=0 );
            assert( idx[i]<lp_cols(lp) );
        }


        int grbError = GRBsetdblattrlist( lp->lp, GRB_DBL_ATTR_START, nNz, idx, coef );
        lp_check_for_grb_error( env, grbError, __FILE__, __LINE__ );
        lp->nModelChanges++;
    }
}
#endif

#ifdef CPX
void lp_config_cpx_params(LinearProgram *lp)
{
    if (lp->maxSeconds != INT_NOT_SET)
        CPXsetdblparam(LPcpxDefaultEnv, CPX_PARAM_TILIM, lp->maxSeconds);
    if (lp->maxSolutions == 1)
        CPXsetdblparam(LPcpxDefaultEnv, CPX_PARAM_EPGAP, 1.0);
    if (lp->maxNodes != INT_NOT_SET)
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_NODELIM, lp->maxNodes);
    if ((lp->silent) || (!lp->printMessages))
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_SCRIND, CPX_OFF);
    else
        CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_SCRIND, CPX_ON);

    if ( lp->mipEmphasis != LP_ME_DEFAULT )
    {
    switch (lp->mipEmphasis)
    {
        case LP_ME_OPTIMALITY:
        {
            CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BESTBOUND );
            break;
        }
        case LP_ME_FEASIBILITY:
        {
            CPXsetintparam(LPcpxDefaultEnv, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_FEASIBILITY );
            break;
        }
    }

    }

    if ( lp->absMIPGap != DBL_MAX )
    {
        double agap;  CPXgetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_AbsMIPGap, &agap );
        printf("changing absolute MIP GAP from %g to %g\n", agap, lp->absMIPGap );
        CPXsetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_AbsMIPGap, lp->absMIPGap );
    }

    if ( lp->relMIPGap != DBL_MAX )
    {
        double rgap;  CPXgetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_MIPGap, &rgap );
        printf("changing relative MIP GAP from %g to %g\n", rgap, lp->relMIPGap );
        CPXsetdblparam( LPcpxDefaultEnv, CPXPARAM_MIP_Tolerances_MIPGap, lp->relMIPGap );
    }

    if ( strlen( lp->solInFN ) )
    {
#define STR_SIZE 256
        vector< pair< string, double > > &mipStart = *(lp->_mipStart);
        char **colNames = (char**) malloc( sizeof(char*)*lp_cols(lp) );
        colNames[0] = (char*) calloc( lp_cols(lp)*STR_SIZE, sizeof(char) );
        for ( int i=1 ; i<lp_cols(lp) ; ++i )
            colNames[i] = colNames[i-1]+STR_SIZE;
        double *colValues = (double*) calloc( lp_cols(lp), sizeof(double) );

        int rCols = lp_read_mip_start( lp->solInFN, colNames, colValues );
        mipStart.clear();
        for ( int i=0 ; i<rCols ; ++i )
            mipStart.push_back( pair< string, double >( colNames[i], colValues[i] ) );

        free( colValues );
        free( colNames[0] );
        free( colNames );
#undef STR_SIZE
        lp->msIdx->clear(); lp->msVal->clear();
    }


    if ( (lp->_mipStart) && (lp->_mipStart->size()) )
        lp_gen_mipstart(lp);

    if ( (lp->msIdx->size()) )
    {
        printf("setting cpx mips start\n");

        int beg[] = { 0 };
        const int effort[] = { CPX_MIPSTART_SOLVEMIP  };

        const int *idx = &((*(lp->msIdx))[0]);
        const double *coef = &((*(lp->msVal))[0]);

        const int nNz = lp->msIdx->size();
        /* checking contents of msIdx */
        for ( int i=0 ; i<nNz ; ++i )
        {
            assert( idx[i]>=0 );
            assert( idx[i]<lp_cols(lp) );
        }

        int cpxError = CPXaddmipstarts( LPcpxDefaultEnv, lp->cpxLP, 1, nNz, beg, idx, coef,
                effort, NULL );
        lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    }
}
#endif // CPX


void lp_write_sol(LinearProgram *lp, const char *fileName)
{
    assert(lp != NULL);

    if (lp_obj_value(lp) == DBL_MAX) {
        fprintf(stderr, "No solution to write.\n");
        abort();
    }

    FILE *fsol = fopen(fileName, "w");
    if (fsol == NULL) {
        fprintf(stderr, "Could not open file %s.\n", fileName);
        abort();
    }

    if (lp->status == LP_OPTIMAL)
        fprintf(fsol, "Optimal (within gap tolerance) - objective value %g - solution time %g\n", lp_obj_value(lp), lp->solutionTime );
    else
        fprintf(fsol, "Stopped on iterations - objective value %g - solution time %g\n", lp_obj_value(lp), lp->solutionTime );

    const int nCols = lp_cols(lp);
    const double *x = lp_x(lp);
    const double *obj = lp_obj_coef( lp );
    char cName[256];
    for (int i = 0 ; (i < nCols) ; ++i) {
        double xv = x[i];

        if (fabs(xv)<EPS)
            continue;

        if (lp_is_integer(lp, i))
            xv = floor(xv + 0.5);

        fprintf(fsol, "%d %s %g %g\n", i, lp_col_name(lp, i, cName), xv, obj[i] );
    }
    fclose(fsol);
}

void lp_parse_options(LinearProgram *lp, int argc, const char **argv)
{
    for (int i = 0 ; (i < argc) ; ++i) {
        if (argv[i][0] != '-')
            continue;

        char optLower[256];
        strncpy(optLower, argv[i], 256);
        int len = strlen(optLower);
        for (int j = 0 ; (j < len) ; ++j)
            optLower[j] = tolower(optLower[j]);

        if (strstr(optLower, "allinteger")) {
            lp->allInteger = 1;
            printf("solving as a pure integer program\n");
            continue;
        }
        if (strstr(optLower, "nomip")) {
            lp->optAsContinuous = 1;
            printf( "solving only root node relaxation.\n" );
            continue;
        }
        if (strstr(optLower, "maxsec")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter number of seconds.\n");
                exit(EXIT_FAILURE);
            }

            int sec = atoi(argv[i + 1]);

            printf("> setting max seconds to %d\n", sec);

            lp_set_max_seconds(lp, sec);
            continue;
        }
        if (strstr(optLower, "mipemphasis")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter MIP emphasis.\n");
                exit(EXIT_FAILURE);
            }

            if ( strcmp(argv[i + 1], "optimality" )==0 )
            {
                lp_set_mip_emphasis( lp, LP_ME_OPTIMALITY );
                printf("MIP emphasis set to optimality.\n");
                continue;
            }
            if ( strcmp(argv[i + 1], "feasibility" )==0 )
            {
                lp_set_mip_emphasis( lp, LP_ME_FEASIBILITY );
                printf("MIP emphasis set to feasibility.\n");
                continue;
            }

            fprintf( stderr, "MIP emphasis not valid: %s. Valid values are: {optimality,feasibility}. ", argv[i+1]);
            exit(EXIT_FAILURE);
        }
        if (strstr(optLower, "maxsol")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter number of solutions.\n");
                exit(EXIT_FAILURE);
            }

            int sol = atoi(argv[i + 1]);

            printf("> setting max solutions to %d\n", sol);

            lp_set_max_solutions(lp, sol);
            continue;
        }
        if (strstr(optLower, "absgap")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter the allowed absolute MIP gap.\n");
                exit(EXIT_FAILURE);
            }

            double agap = atof(argv[i + 1]);

            lp_set_abs_mip_gap( lp, agap );
            continue;
        }
        if (strstr(optLower, "relgap")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter the relative MIP gap.\n");
                exit(EXIT_FAILURE);
            }

            double rgap = atof(argv[i + 1]);

            lp_set_rel_mip_gap( lp, rgap );
            continue;
        }
        if (strstr(optLower, "soloutfn")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter solution file name.\n");
                exit(EXIT_FAILURE);
            }

            printf("> setting solution output file name to %s\n", argv[i + 1]);

            lp_set_sol_out_file_name(lp, argv[i + 1]);
            continue;
        }
        if (strstr(optLower, "solinfn")) {
            if (i + 1 == argc) {
                fprintf(stderr, "enter solution file name.\n");
                exit(EXIT_FAILURE);
            }

            printf("> setting solution input file name to %s\n", argv[i + 1]);

            lp_set_sol_in_file_name(lp, argv[i + 1]);
            continue;
        }
    }
}

void lp_set_sol_out_file_name(LinearProgram *lp, const char *sfn)
{
    assert(lp);

    strcpy(lp->solOutFN, sfn);
}

void lp_set_sol_in_file_name(LinearProgram *lp, const char *sfn)
{
    assert(lp);

    strcpy(lp->solInFN, sfn);
}

void lp_load_mip_starti( LinearProgram *lp, int count, const int *colIndexes, const double *colValues )
{
    char colName[256];
    lp->_mipStart->clear();
    for (int i = 0 ; (i < count) ; i++)
        lp->_mipStart->push_back(pair< string, double >( lp_col_name( lp,colIndexes[i],colName ), colValues[i]) );
}

void lp_load_mip_start(LinearProgram *lp, int count, char **colNames, const double *colValues)
{
    lp->_mipStart->clear();
    for (int i = 0 ; (i < count) ; i++)
        lp->_mipStart->push_back(pair< string, double >(colNames[i], colValues[i]));
}

void lp_chg_obj(LinearProgram *lp, int count, int idx[], double obj[])
{
#ifdef GRB
    int grbError = GRBsetdblattrlist( lp->lp, "Obj", count, idx, obj );

    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif
#ifdef CPX
    int cpxError = CPXchgobj(LPcpxDefaultEnv, lp->cpxLP, count, idx, obj);

    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    for (int i=0 ; (i<count) ; i++ )
        lp->_lp->setObjCoeff( idx[i], obj[i]);
#endif
}

int lp_nz(LinearProgram *lp)
{
    assert(lp);
#ifdef CPX
    return CPXgetnumnz(LPcpxDefaultEnv, lp->cpxLP);
#endif
#ifdef CBC
    return lp->_lp->getNumElements();
#endif
#ifdef GLPK
    return glp_get_num_nz( lp->_lp );
#endif
#ifdef GRB
    int nzs = 0;
    int grbError = GRBgetintattr( lp->lp, "NumNZs", &nzs );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    return nzs;
#endif
    return 0;
}

void lp_cols_by_type( LinearProgram *lp, int *binaries, int *integers, int *continuous )
{
    *binaries = *integers = *continuous = 0;

    for (int i=0 ; (i<lp_cols(lp)) ; i++ )
    {
        const double lb = lp_col_lb( lp, i );
        const double ub = lp_col_ub( lp, i );

        if (lp_is_integer( lp, i ))
        {
            if ( ( fabs( lb ) <= EPS ) && (fabs( ub - 1.0 ) <= EPS) )
                (*binaries)++;
            else
                (*integers)++;
        }
        else
            (*continuous)++;
    }
}

void lp_help_options()
{
    printf("options:\n");
    printf("\t-nomip             :  solves only the linear programming relaxation\n");
    printf("\t-allinteger        :  solves as a pure integer program\n");
    printf("\t-maxSec sec        :  specifies timelimit of 'sec' seconds\n");
    printf("\t-maxSol sol        :  search will be ended when 'sol' solutions are found.\n");
    printf("\t-absgap gap        :  set allowed the absolute MIP gap to 'gap'.\n");
    printf("\t-relgap gap        :  set allowed the relative MIP gap to 'gap'.\n");
    printf("\t-solOutFN solfname :  file name to save solution\n");
    printf("\t-solInFN  solfname :  file name to read solution\n");
}

int lp_read_mip_start( const char *fileName, char **colNames, double colValues[] )
{
#define STR_SIZE 256
    FILE *f = fopen( fileName, "r" );
    if (!f)
    {
        fprintf( stderr, "Could not open mipstart from file %s at %s:%d.\n", fileName, __FILE__, __LINE__ );
        abort();
    }
    char line[STR_SIZE];

    int nLine = 0;
    int nRCols = 0;
    while (fgets( line, STR_SIZE, f ))
    {
        ++nLine;
        char col[4][STR_SIZE];
        int nread = sscanf( line, "%s %s %s %s", col[0], col[1], col[2], col[3] );
        if (!nread)
            continue;
        /* line with variable value */
        if (strlen(col[0])&&isdigit(col[0][0])&&(nread>=3))
        {
            if (!lp_str_is_num(col[0]))
            {
                fprintf( stderr, "LP warning: reading: %s, line %d - first column in mipstart file should be numeric, ignoring.", fileName, nLine );
                continue;
            }
            if (!lp_str_is_num(col[2]))
            {
                fprintf( stderr, "LP warning: reading: %s, line %d - third column in mipstart file should be numeric, ignoring.", fileName, nLine  );
                continue;
            }

            //int idx = atoi( col[0] );
            char *name = col[1];
            double value = atof( col[2] );
            //double obj = 0.0;
            //         if (nread >= 4)
            //            obj = atof( col[3] );

            strncpy( colNames[nRCols], name, STR_SIZE );
            colValues[nRCols] = value;

            ++nRCols;
        }
    }

    if (nRCols)
        printf("LP: mipstart values read for %d variables.\n", nRCols );
    else
        printf("LP: no mipstart solution read from %s.\n", fileName );

    fclose(f);

    return nRCols;
#undef STR_SIZE
}

void lp_set_abs_mip_gap( LinearProgram *lp, const double _value )
{
    lp->absMIPGap = _value;
}

void lp_set_rel_mip_gap( LinearProgram *lp, const double _value )
{
    lp->relMIPGap = _value;
}

bool lp_str_is_num( const char *str )
{
    const size_t l = strlen(str);

    for ( size_t i=0 ; i<l ; ++i )
        if (!(isdigit(str[i])||(str[i]=='.')))
            return false;

    return true;
}

char lp_is_binary( LinearProgram *lp, const int j )
{
    LP_CHECK_COL_INDEX( lp, j );

    return ( (lp_is_integer(lp,j)) && (fabs(lp_col_lb(lp,j))<= EPS)  && (fabs(lp_col_ub(lp,j)-1.0)<=EPS) );
}

int lp_row_type( LinearProgram *lp, const int row )
{
    LP_CHECK_ROW_INDEX( lp, row );

    lp->_idx->resize( lp_cols(lp) );
    lp->_coef->resize( lp_cols(lp) );

    int *idx = &((*lp->_idx)[0]);
    double *coef = &((*lp->_coef)[0]);

    int result = CONS_OTHER;

    double minc = DBL_MAX, maxc = -DBL_MAX;
    int nz = lp_row( lp, row, idx, coef );
    assert( nz>=0 );

    int nbinaries = 0;
    int nint = 0;
    int ncont = 0;


    int nIntCoef = 0;
    int j;
    for ( j=0 ; (j<nz) ; ++j )
    {
        minc = std::min( minc, coef[j] );
        maxc = std::max( maxc, coef[j] );
        if (lp_is_binary(lp,idx[j]))
            nbinaries++;
        else
            if (lp_is_integer(lp,idx[j]))
                nint++;
            else
                ncont++;
        nIntCoef += is_integer( coef[j] );
    }
    char sense = lp_sense( lp, row );
    double rhs = lp_rhs( lp, row );
    char allOneLHS = ( (fabs(minc-maxc)<=EPS) && (fabs(minc-1.0)<=EPS) );
    char rhsOne = fabs(rhs-1.0) <= EPS;

    /* binaries, all positive and integral */
    if ( (nbinaries==nz)&&(minc>=0.98)&&(rhs>=0.98) )
    {
        switch (sense)
        {
            case 'E':
            {
                if (allOneLHS )
                {
                    if ( rhsOne )
                        result = CONS_PARTITIONING;
                    else
                        if ( rhs >= 1.99 )
                            result = CONS_CARDINALITY;
                }
                goto RETURN_POINT;
            }
            case 'L':
            {
                if ( (allOneLHS) && (rhsOne) )
                {
                    result = CONS_PACKING;
                    goto RETURN_POINT;
                }
                else
                if (allOneLHS)
                {
                    result = CONS_INV_KNAPSACK;
                    goto RETURN_POINT;
                }
                else
                if ( (maxc>=1.99) && (rhs>=1.99) )
                {
                    result = CONS_KNAPSACK;
                    goto RETURN_POINT;
                }

               goto RETURN_POINT;
            }
            case 'G':
            {
                if( (allOneLHS) && (rhsOne) )
                {
                    result = CONS_COVERING;
                    goto RETURN_POINT;
                }
                goto RETURN_POINT;

            }
        }
    }

    /*  = 0,  flow constraints */
    if ((fabs(rhs)<=1e-8) && (sense=='E'))
    {
        if ( ( minc <= -0.98 ) && (maxc >= 0.98 ) )
        {
            if ( nbinaries == nz )
            {
                result = CONS_FLOW_BIN;
                goto RETURN_POINT;
            }
            else
            {
                if ( nint+nbinaries == nz )
                {
                    result = CONS_FLOW_INT;
                    goto RETURN_POINT;
                }
                else
                {
                    result = CONS_FLOW_MX;
                    goto RETURN_POINT;
                }
            }
        }
    }

RETURN_POINT:

    return result;
}

double lp_best_bound( LinearProgram *lp )
{
    return lp->bestBound;
}

void lp_rows_by_type( LinearProgram *lp, int rtype[] )
{
    memset( rtype, 0, sizeof(int)*CONS_NUMBER );
    int i;
    for ( i=0 ; (i<lp_rows(lp)) ; ++i )
        rtype[lp_row_type(lp,i)]++;
}

void lp_set_integer( LinearProgram *lp, int nCols, int cols[] )
{
#ifdef CBC
    printf("transforming %d variables into integer ones.\n", nCols );
    lp->osiLP->setInteger( cols, nCols );
#endif
#ifdef GLPK
    for ( int i=0 ; i<nCols ; ++i )
    {
    if ((lp_col_lb(lp,i)<=-1e-5)||(lp_col_ub(lp,i)>=1.0+1e-5))
            glp_set_col_kind( lp->_lp, cols[i]+1, GLP_IV);
        else
            glp_set_col_kind( lp->_lp, cols[i]+1, GLP_BV);
    }
#endif
#ifdef GRB
    char *vType = new char[nCols];

    memset( vType, GRB_INTEGER, sizeof(char)*nCols );
    int grbError = GRBsetcharattrlist( lp->lp, "VType", nCols, cols, vType );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );

    delete[] vType;
#endif
}

char *lp_status_str( int status, char *statusStr )
{
    switch (status)
    {
    case LP_OPTIMAL:
        sprintf( statusStr, "LP_OPTIMAL");
        break;
    case LP_UNBOUNDED:
        sprintf( statusStr, "LP_UNBOUNDED");
        break;
    case LP_INFEASIBLE:
        sprintf( statusStr, "LP_INFEASIBLE");
        break;
    case LP_FEASIBLE:
        sprintf( statusStr, "LP_FEASIBLE");
        break;
    case LP_INTINFEASIBLE:
        sprintf( statusStr, "LP_INTINFEASIBLE");
        break;
    case LP_NO_SOL_FOUND:
        sprintf( statusStr, "LP_NO_SOL_FOUND");
        break;
    case LP_ERROR:
        sprintf( statusStr, "LP_ERROR");
        break;
    default:
        fprintf( stderr, "lp status not recognized: %d\n", status );
        abort();
    }

    return statusStr;
}

int *lp_original_colummns( LinearProgram *lp )
{
    return &((*lp->_orig)[0]);
#ifdef GRB
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif
}

void lp_mipstart_debug( LinearProgram *lp )
{
    vector< double > x( lp_cols(lp), 0.0 );
#ifdef GRB
    if ( lp->tmpCols || lp->tmpRows || lp->nModelChanges )
    {
        int grbError = GRBupdatemodel( lp->lp );
        lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
        lp->tmpCols = lp->tmpRows = lp->nModelChanges = 0;
    }

#endif

    vector< pair< string, double > > &mipStart = *(lp->_mipStart);
    if ( strlen( lp->solInFN ) )
    {
#define STR_SIZE 256
        char **colNames = (char**) malloc( sizeof(char*)*lp_cols(lp) );
        colNames[0] = (char*) calloc( lp_cols(lp)*STR_SIZE, sizeof(char) );
        for ( int i=1 ; i<lp_cols(lp) ; ++i )
            colNames[i] = colNames[i-1]+STR_SIZE;
        double *colValues = (double*) calloc( lp_cols(lp), sizeof(double) );

        int rCols = lp_read_mip_start( lp->solInFN, colNames, colValues );
        mipStart.clear();
        for ( int i=0 ; i<rCols ; ++i )
            mipStart.push_back( pair< string, double >( colNames[i], colValues[i] ) );

        free( colValues );
        free( colNames[0] );
        free( colNames );
#undef STR_SIZE
        lp->msIdx->clear(); lp->msVal->clear();
    }

    if ( (lp->_mipStart) && (lp->_mipStart->size()) )
    {
        lp_gen_mipstart(lp);

        const vector< int > &msIdx = *(lp->msIdx);
        const vector< double > &msVal = *(lp->msVal);

        assert( msIdx.size() == msVal.size() );
        int nFix = (int) msIdx.size();
        for ( int i=0 ; (i<nFix) ; ++i )
        {
            int idxX = msIdx[i];
            assert( msIdx[i]>=0 && msIdx[i]<(int)x.size() );
            x[ idxX ] = msVal[ i ];
            double lbCol = lp_col_lb( lp, idxX );
            double ubCol = lp_col_ub( lp, idxX );
            if ( msVal[i]<lbCol-1e-10 || msVal[i]>ubCol+1e+10 )
            {
                char colName[256];
                lp_col_name( lp, idxX, colName );
                fprintf( stderr, "variable %d (%s) has bound [%g,%g] but value in mipstart is %g\n", idxX, colName, lbCol, ubCol, x[idxX]  );
                abort();
            }
        }
    }
    else
    {
        printf("No MIPStart informed for debugging.\n");
        return;
    }

    vector< int > idx( lp_cols(lp), 0 ) ;
    vector< double > coef( lp_cols(lp), 0 ) ;

    vector< char > ivActive( lp_cols( lp ), 0 );
    vector< int >::iterator msIt = lp->msIdx->begin();
    for ( ; (msIt != lp->msIdx->end()) ; ++msIt )
        if (lp_is_binary( lp, *msIt ))
            ivActive[*msIt]++;



    bool bugged = false;

    for ( int i=0 ; i<lp_rows(lp) ; ++i )
    {
        double rhs = lp_rhs( lp, i );
        const char sense = lp_sense( lp, i );
        int nzRow = lp_row( lp, i, &(idx[0]), &(coef[0]) );
        double minRow = DBL_MAX;
        double maxRow = -DBL_MAX;
        if ( dbl_equal( rhs, 1.0 ) )
        {
            for ( int j=0 ; (j<nzRow) ; ++j )
            {
                minRow = std::min( minRow, coef[j] );
                maxRow = std::max( maxRow, coef[j] );
            }

            if ( dbl_equal( minRow, 1.0 ) && dbl_equal( maxRow, 1.0) )
            {
                /*
                {
                    char rowN[256];
                    printf("row %s\n", lp_row_name( lp, i, rowN ));
                }*/

                int nActive = 0;
                for ( int k=0 ; k<nzRow ; ++k )
                    nActive += ivActive[idx[k]];


                if ( nActive>=2 && ( sense=='L' || sense == 'E' ) )
                {
                    int nPrints = 0;
                    char confNames[2][256] = { "","" };
                    for ( int j=0 ; ((j<nzRow)&&(nPrints<2)) ; ++j )
                        if (ivActive[idx[j]])
                        {
                            char colName[256];
                            strcpy( confNames[nPrints], lp_col_name( lp, idx[j], colName ) );
                            nPrints++;
                        }

                    char rowName[256];
                    fprintf( stderr, "MIPStart: packing/partitioning constraint %s has at least two active variables: %s and %s\n", lp_row_name(lp, i, rowName), confNames[0], confNames[1] );
                    bugged = true;
                }
            } // minRow e maxRow
            else
            {
                // positive coefficients >= positive
                if ( minRow>0 && rhs>0 && sense=='G' )
                {
                    double lhs = 0;
                    for ( int k=0 ; k<nzRow ; ++k )
                        lhs += coef[ k ] * x[ idx[k] ];

                    if ( lhs <= rhs-1e-10 )
                    {
                        char rowName[256];
                        lp_row_name( lp, i, rowName );
                        fprintf( stderr, "Constraint %s (lhs >= %g) has lhs from mipstart %g.\n", rowName, rhs, lhs );
                        bugged = true;
                    }
                }
            } // not packing or partitioning
        } // rhs 1.0
    }

    if (bugged)
        abort();
}

void lp_remove_row( LinearProgram *lp, int idxRow )
{
    LP_CHECK_ROW_INDEX( lp, idxRow );
#ifdef GRB
    int pIdxRow = idxRow;
    int grbError = GRBdelconstrs( lp->lp, 1, &pIdxRow );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
    lp->nModelChanges++;
#endif
#ifdef CPX
    int cpxError = CPXdelrows( LPcpxDefaultEnv, lp->cpxLP, idxRow, idxRow );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
#endif
#ifdef CBC
    int idx = idxRow;
    lp->osiLP->deleteRows(1, &idx );
#endif
#ifdef GLPK
    int idx = idxRow+1;
    glp_del_rows( lp->_lp, 1, &idx );
#endif
}

bool dbl_equal( const double v1, const double v2 )
{
    return (fabs( v2-v1 ) <= 1e-9);
}

double round( const double v )
{
   return (double) ( (long long) (v+0.5) );
}

bool is_integer( const double v )
{
    if ( (v==DBL_MAX) || (v==DBL_MIN) )
        return true;

    return (  fabs(round(v) - v) <= 1e-10 );
}

void lp_close_env()
{
#ifdef GRB
    if (LPgrbDefaultEnv != NULL)
    {
        GRBfreeenv( LPgrbDefaultEnv );
    }
#endif
#ifdef CPX
    if ( LPcpxDefaultEnv != NULL )
    {
        CPXcloseCPLEX( &LPcpxDefaultEnv );
    }
#endif
}

double* lp_reduced_cost(LinearProgram* lp)
{
    assert(lp != NULL);

    if (lp->nOptimizations == 0) {
        fprintf(stderr, "No optimizations have been made with this model.\n");
        abort();
    }

    if (lp->status != LP_OPTIMAL) {
        fprintf(stderr, "\n\nERROR: no dual solution available.\n At: %s:%d\n\n", __FILE__, __LINE__);
        abort();
    }

    return &((*(lp->_rc))[0]);
}

void lp_set_store_names(int store)
{
    LPstoreNames = store;
}


/* Define my callback function */
int getcallback(LinearProgramPtr lp, int time){
    int erro = GRBsetcallbackfunc(lp->lp, mycallback, (void*) &time);

}

int setNodeInfinity(LinearProgramPtr lp){
#ifdef GRB
        GRBenv *env = GRBgetenv( lp->lp );
        GRBsetdblparam(env, GRB_DBL_PAR_NODELIMIT, GRB_INFINITY );
#endif // GRB
}

int mycallback(GRBmodel *model, void *cbdata, int where, void *usrdata)
{
#ifdef GRB
    GRBenv *modelenv = GRBgetenv(model);
    int time;
    double runTime;

    time = *((int*)usrdata);
    if (where == GRB_CB_MIPSOL) {
        double node;
        int cont;
        GRBcbget(cbdata, where, GRB_CB_MIPSOL_NODCNT, &node);
        GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOLCNT, &cont);
        GRBcbget(cbdata, where, GRB_CB_RUNTIME, &runTime);
        time = time - (int)runTime;
        printf("%d\n", time);
        if(cont>1){
            GRBsetdblparam(modelenv, GRB_DBL_PAR_NODELIMIT, node+500 );
            GRBsetdblparam( modelenv, GRB_DBL_PAR_TIMELIMIT, time);
        }
    }
    return 0;
#endif // GRB
}


void lp_save_mip_start( LinearProgram *lp, const char *fileName )
{
    char fname[512];
    strcpy( fname, fileName );
    if (strstr(fname, ".")==NULL)
        strcat(fname, ".sol");

    assert( lp );
    assert( lp->_mipStart );

    const std::vector< std::pair< std::string, double > > &mipStart = *lp->_mipStart;
    assert( mipStart.size() );

    FILE *f = fopen( fname, "w" );
    fprintf( f, "Feasible - 0.00\n" );
    size_t nv = mipStart.size(), i=0;
    for ( i=0 ; (i<nv) ; ++i )
    {
        const string &str = mipStart[i].first;
        const double &val = mipStart[i].second;
        fprintf( f, "%zu %s %g\n", i, str.c_str(), val );
    }
    fclose(f);
}

void lp_remove_rows( LinearProgram *lp, int nRows, int *rows )
{
    assert( lp );
    for ( int i=0 ; (i<nRows); ++i )
    {
        LP_CHECK_ROW_INDEX( lp, rows[i] );
    }
// more expensive check
#ifdef DEBUG
    for ( int i1=0 ; (i1<nRows); ++i1 )
    {
        for ( int i2=i1+1 ; (i2<nRows); ++i2 )
        {
            // no repeated entries
            assert( rows[i1]!=rows[i2] );
        }
    }
#endif
#ifdef CBC
    lp->osiLP->deleteRows( nRows, rows );
#endif
#ifdef GRB
    int grbError =
        GRBdelconstrs( lp->lp, nRows, rows );
    lp_check_for_grb_error( LPgrbDefaultEnv, grbError, __FILE__, __LINE__ );
#endif
#ifdef CPX
    int *rset = (int*) calloc( sizeof(int), lp_rows(lp) );
    for ( int i=0 ; (i<nRows) ; ++i )
        rset[rows[i]] = 1;
    int cpxError =
        CPXdelsetrows( LPcpxDefaultEnv, lp->cpxLP, rset );
    lp_check_for_cpx_error(LPcpxDefaultEnv, cpxError, __FILE__, __LINE__);
    free( rset );
#endif
#ifdef GLPK
    fprintf( stderr, "Call not implemented in LP yet.\n");
    abort();
#endif
}

