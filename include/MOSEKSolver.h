

#ifndef MOSEKSOLVER_H
#define MOSEKSOLVER_H


#include "mosek.h"
#include "LPSolver.h"
#include <vector>
#include <iostream>

static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */



class MOSEKSolver : public LPSolver{

  private:
  
    typedef LPSolver::Status Status;

    MSKenv_t     env;
    MSKtask_t    task;
    double *primal;
    double *dual;
    MSKint32t optimType;
    std::vector< MSKstakeye > colStatus;
    std::vector< MSKstakeye > rowStatus;
    MSKrescodee res;


  public:

    MOSEKSolver( MSKint32t optimizer = MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX ) {
        MSK_makeenv(&env,NULL);
        task = NULL;
        primal = NULL;
        dual=NULL;
        optimType=optimizer;
    };

    ~MOSEKSolver(){
      deleteLP();
      MSK_deleteenv(&env);
    };


   virtual void solveLP(){

     if(primal != NULL){
       delete[] primal;
     }
     if(dual != NULL){
       delete[] dual;
     }
     primal = new double[ getNumberOfColumns() ];
     dual = new double[ getNumberOfRows() ];


     MSK_putskx (
         task,
         MSK_SOL_BAS,
         colStatus.data() );

     MSK_putskc (
         task,
         MSK_SOL_BAS,
         rowStatus.data() );


     res = MSK_optimize(task);
     std::cout << "MSK result code: " << res << std::endl;

     MSK_getxx(task,
         MSK_SOL_BAS,    /* Request the basic solution. */
         primal);

     MSK_gety(task,
         MSK_SOL_BAS,    /* Request the basic solution. */
         dual);

     MSK_getskx(task, MSK_SOL_BAS, colStatus.data() );
     MSK_getskc(task, MSK_SOL_BAS, rowStatus.data() );
   };

      virtual double getObjectiveValue(){
     double obj;
     MSK_getprimalobj(task, MSK_SOL_BAS, &obj);
     return obj;
   };


   virtual long getIterationCount(){
     static long n=16;
     static MSKiinfiteme iters[]= {
       MSK_IINF_INTPNT_ITER,
       MSK_IINF_SIM_DUAL_DEG_ITER,
       MSK_IINF_SIM_DUAL_INF_ITER,
       MSK_IINF_SIM_DUAL_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_DEG_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_INF_ITER,
       MSK_IINF_SIM_NETWORK_DUAL_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_DEG_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_INF_ITER,
       MSK_IINF_SIM_NETWORK_PRIMAL_ITER,
       MSK_IINF_SIM_PRIMAL_DEG_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_DEG_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_INF_ITER,
       MSK_IINF_SIM_PRIMAL_DUAL_ITER,
       MSK_IINF_SIM_PRIMAL_INF_ITER,
       MSK_IINF_SIM_PRIMAL_ITER
     };

     long nIter = 0;

     for(long i=0; i<n; i++){
       MSKint32t count;

       MSK_getintinf (
           task,
           iters[i],
           &count
           );
       nIter += count;
     }
     return nIter;

   };


   virtual long getNumberOfRows(){
     MSKint32t numcon;
     MSK_getnumcon (
         task,
         &numcon);
     return numcon;
   };

   virtual long getNumberOfColumns(){
     MSKint32t numvar;
     MSK_getnumvar (
         task,
         &numvar);
     return numvar;
   };

   virtual void setupStandardBasis(){
     for(long i= 0; i< getNumberOfColumns(); i++){
       colStatus[i] = MSK_SK_LOW;
     }
     for(long i= 0; i< getNumberOfRows(); i++){
       rowStatus[i] = MSK_SK_BAS;
     }
   };



   virtual void createLP(long nSource, long nTarget){
     deleteLP();

     MSK_maketask(env, 0, 0,&task);
     MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
     MSK_putintparam(task, MSK_IPAR_OPTIMIZER,
     //    MSK_OPTIMIZER_FREE_SIMPLEX);
        optimType);
     MSK_putintparam(task, MSK_IPAR_SIM_HOTSTART, MSK_SIM_HOTSTART_STATUS_KEYS);
     MSK_putintparam(task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
     //MSK_putintparam(task, MSK_IPAR_PRESOLVE_ELIMINATOR_USE, MSK_PRESOLVE_MODE_OFF);
     MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL, printstr);
   };



   virtual void addColumns(long n){
     MSK_appendvars(task,n);
     for(long i=0; i<n; i++){
       colStatus.push_back( MSK_SK_LOW );
     }
   };



   virtual void addRows(long n){
     MSK_appendcons(task, n);
     for(long i=0; i<n; i++){
       rowStatus.push_back( MSK_SK_LOW );
     }
   };



   virtual double getRowDual(long row){
     return dual[row];
   };



   virtual double getColumnPrimal(long col){
     return primal[col];
   };



   virtual void setRowBounds(long i, double mass){

     MSK_putconbound(task,
                          i,           /* Index of constraint.*/
                          MSK_BK_FX,      /* Bound key.*/
                          mass,      /* Numerical value of lower bound.*/
                          mass);     /* Numerical value of upper bound.*/
   };

   virtual double getRowBounds(long i){
     MSKboundkeye bk;
     MSKrealt bl;
     MSKrealt bu;
     MSK_getbound (
         task,
         MSK_ACC_CON,
         i,
         &bk,
         &bl,
         &bu);
     return bu;

   };

   virtual void setColumnBounds(long i, double lb, double ub){
     MSK_putvarbound(task,
         i,           /* Index of variable.*/
         MSK_BK_RA,      /* Bound key.*/
         lb,      /* Numerical value of lower bound.*/
         ub);     /* Numerical value of upper bound.*/
   };

   virtual void setColumnBoundsLower(long i, double lb){
     MSK_putvarbound(task,
         i,           /* Index of variable.*/
         MSK_BK_LO,      /* Bound key.*/
         lb,      /* Numerical value of lower bound.*/
         MSK_SK_INF);     /* Numerical value of upper bound.*/
   };

   virtual void setColumnObjective(long i, double cost){
     MSK_putcj(task, i, cost);
   };

   virtual void setColumnCoefficients( long col, long s, long t){
     MSKint32t ind[2] = {(int)s,(int)t};
     double val[2] = {1,-1};
     MSK_putacol(task,
                        col,                 /* Variable (column) index.*/
                        2, /* Number of non-zeros in column j.*/
                        ind,     /* Pointer to row indexes of column j.*/
                        val);    /* Pointer to Values of column j.*/
   };

   virtual long getColumn(long col, long *ind, double *val){
     MSKint32t colid = (int) col;
     MSKint32t n = 0;
     MSKint32t indi[2];
     MSK_getacol(task, colid, &n, indi, val);
     ind[0] = indi[0];
     ind[1] = indi[1];
     return n;
   };

   virtual Status getColumnStatus(long col){
     return convertFromMOSEK( colStatus[col] );
   };

   virtual Status getRowStatus(long row){
     return convertFromMOSEK( rowStatus[row] );
   };

   virtual void setColumnStatus(long col, Status s){
     colStatus[col] = convertToMOSEK(s);
   };

   virtual void setRowStatus(long row, Status s){
     rowStatus[row] = convertToMOSEK(s);
   };


   virtual bool isOptimal(){
     return res == MSK_RES_OK;
   };

  private:

   
   virtual void deleteLP(){
     if(task != NULL){
       MSK_deletetask(&task);
       task = NULL;
       delete[] primal;
       delete[] dual;
       colStatus.clear();
       rowStatus.clear();
       dual = NULL;
       primal = NULL;
     }
   };




   Status convertFromMOSEK(MSKstakeye s){
     switch(s){
       case MSK_SK_BAS:
         return LPSolver::BASIC;
       case MSK_SK_UPR:
         return LPSolver::UPPER;
       case MSK_SK_LOW:
         return LPSolver::LOWER;
       case MSK_SK_UNK:
         return LPSolver::UNKNOWN;
       case MSK_SK_SUPBAS:
         return LPSolver::SUPERBASIC;
       case MSK_SK_FIX:
         return LPSolver::FIXED;
       case MSK_SK_INF:
         return LPSolver::INF;
       case MSK_SK_END:
         return LPSolver::END;
     }
     return LPSolver::END;
   };


   MSKstakeye convertToMOSEK(Status s){
     switch(s){
       case LPSolver::BASIC:
         return MSK_SK_BAS;
       case LPSolver::UPPER:
         return MSK_SK_UPR;
       case LPSolver::LOWER:
         return MSK_SK_LOW;
       case LPSolver::UNKNOWN:
       case LPSolver::FREE:
         return MSK_SK_UNK;
       case LPSolver::SUPERBASIC:
         return MSK_SK_SUPBAS;
       case LPSolver::FIXED:
         return MSK_SK_FIX;
       case LPSolver::INF:
         return MSK_SK_INF;
       case LPSolver::END:
         return MSK_SK_END;

     }
     return MSK_SK_END;
   };



};


#endif
