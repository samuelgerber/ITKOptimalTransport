#ifndef LPSOLVER_H
#define LPSOLVER_H

#include "MultiscaleTransport.h"

class LPSolver {
  public:


    enum Status { UPPER, LOWER, BASIC, UNKNOWN, 
                  SUPERBASIC, FIXED, INF, END, FREE};




  public:



    virtual ~LPSolver(){};

    virtual void solveLP() = 0;
    virtual bool isOptimal() = 0;
    virtual long getNumberOfRows() = 0;
    virtual long getNumberOfColumns() = 0;
    virtual double getObjectiveValue() = 0;
    virtual long getIterationCount() = 0;
  
    virtual void setupStandardBasis() = 0;
   
    virtual void createLP(long nSource, long nTarget ) = 0;


    virtual void addColumns(long n) = 0;
    virtual void addRows(long n) = 0;
   
    virtual double getRowDual(long row) = 0;
    virtual double getColumnPrimal(long col) = 0;
   
    virtual void setRowBounds(long row, double mass) = 0;
    virtual double getRowBounds(long row) = 0;
  
    virtual void setColumnBoundsLower(long col, double lb) = 0; 
    virtual void setColumnBounds(long col, double lb, double ub) = 0;
    virtual void setColumnObjective(long col, double cost) = 0;
    virtual void setColumnCoefficients( long col, long s, long t) = 0; 
    
    virtual long getColumn(long col, long *ind, double *val) = 0;
    
    virtual Status getRowStatus(long row) = 0;
    virtual void setRowStatus(long row, Status s) = 0;


    virtual Status getColumnStatus(long col) = 0;
    virtual void setColumnStatus(long col, Status s) = 0;
   



};


#endif


