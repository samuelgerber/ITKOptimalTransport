

#ifndef LEMONSOLVER_H
#define LEMONSOLVER_H

#include "LPSolver.h"

#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>

class LemonSolver : public LPSolver{

  private:
 
    typedef typename LPSolver::Status Status;

    std::vector<long> sInd;
    std::vector<long> tInd;

    std::vector<double> coeff;
    std::vector<double> mass;
    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<double> colLB;
    std::vector<double> colUB;

    std::vector<Status> colStatus;
    std::vector<Status> rowStatus;

    long solstat;
    long iCount;
    double objValue;
    long ns;
    long nt;
    bool success;




  public:

    LemonSolver() {
      ns = 0;
      nt = 0;
      success = false;
    };

    ~LemonSolver(){
      deleteLP();
    };




   virtual void solveLP(){

     using namespace lemon;

     SmartDigraph graph;
     typedef SmartDigraph::Node Node;
     typedef SmartDigraph::NodeIt NodeIt;
     typedef SmartDigraph::Arc Arc;
     typedef SmartDigraph::ArcIt ArcIt;
     typedef SmartDigraph::InArcIt InArcIt;
     typedef SmartDigraph::OutArcIt OutArcIt;
     typedef SmartDigraph::NodeMap<long long> IntNodeMap;
     typedef SmartDigraph::ArcMap<long long> IntArcMap;

     long status;

     //Lemon allwos integer only:
     //Scale cost and capacities
     static const long long maxVal = 10000000000L;

     double costScaling = 0;
     for(long i=0; i<coeff.size(); i++){
       if(coeff[i] > costScaling){
         costScaling = coeff[i];
       }
     }
     if( costScaling == 0){
       costScaling = maxVal;
     }
     else{
       costScaling = ((double) maxVal ) / costScaling;
     }


     double mPositive = 0;
     double mNegative = 0;
     for(long i=0; i < mass.size(); i++){
        double m = mass[i];
        if( m > 0 ){
          mPositive += m;
        }
        else{
          mNegative += m;
        }
     }
     double maxCapacity = std::max( mPositive, -mNegative);
     double capacityScaling = ((double) maxVal ) / maxCapacity;

#ifdef VERBOSE
     std::cout << "capacityScaling: "<< capacityScaling << std::endl;
     std::cout << "costScaling: "<< costScaling << std::endl;
#endif

     //Setup nodes
     long long massPositive = 0;
     long long massNegative = 0;
     graph.reserveNode( mass.size() +1 );
     for(long i=0; i <= mass.size(); i++){
        graph.addNode();
     }
     IntNodeMap supply(graph);

     long long maxMass = 0;
     int maxMassID = -1;
     long long minMass = 0;
     int minMassID = -1;
     for(long i=0; i<mass.size(); i++){
        long long m = (long long) ( mass[i] * capacityScaling );
        //double m = mass[i];
        supply[ graph.nodeFromId(i) ] = m;

        if( m > 0 ){
          massPositive += m;
        }
        else{
          massNegative += m;
        }

        if( m > maxMass ){
          maxMass = m;
          maxMassID = i;
        }
        if( m < minMass ){
          minMass = m;
          minMassID = i;
        }

     }

     long long massImbalance = massPositive + massNegative;
     if(massImbalance > 0 ){
       supply[ graph.nodeFromId(maxMassID) ] = maxMass - massImbalance;
     }
     if(massImbalance < 0 ){
       supply[ graph.nodeFromId(minMassID) ] = minMass - massImbalance;
     }
     
#ifdef VERBOSE
       std::cout << "Mass positive: " << massPositive << std::endl;
       std::cout << "Mass negative: " << massNegative << std::endl;
       std::cout << "Mass imbalance: " << massImbalance << std::endl;
       std::cout << "Max  Mass: " << maxMass<< std::endl;
       std::cout << "Min  Mass: " << minMass<< std::endl;
       std::cout << "Max  Mass ID: " << maxMassID<< std::endl;
       std::cout << "Min  Mass ID: " << minMassID<< std::endl;
#endif

     


     graph.reserveArc( sInd.size()  );


     //Add regular node
     for( long i=0; i < sInd.size(); i++){
       graph.addArc( graph.nodeFromId(sInd[i]), graph.nodeFromId(tInd[i]) );
     }

     IntArcMap capacity(graph), lower(graph), cost(graph);
     for( long i=0; i < coeff.size(); i++){
       Arc a = graph.arcFromId(i);
       if( colUB[i] > maxCapacity ){
         capacity[a] = (long long) (capacityScaling * maxCapacity ) + 1;
       }
       else{
         capacity[a] = (long long) (capacityScaling * colUB[i] ) + 1;
       }
       //capacity[a] = colUB[i];
       lower[a] = (long long) std::max( 0LL, (long long) (capacityScaling * colLB[i] ) -1 );
       //lower[a] = colLB[i];
       cost[a] = (long long)( costScaling * coeff[i] );
       //cost[a] = coeff[i];
     }


     //Solve
     NetworkSimplex<SmartDigraph, long long> simplex(graph);
     simplex.upperMap(capacity);
     simplex.lowerMap(lower);
     simplex.costMap(cost);
     simplex.supplyMap(supply);


     NetworkSimplex<SmartDigraph, long long>::ProblemType res = simplex.run();
     success = res == NetworkSimplex<SmartDigraph, long long>::OPTIMAL;
     objValue  = simplex.totalCost<double>();
     objValue /= capacityScaling;
     objValue /= costScaling;
     iCount = 1;

#ifdef VERBOSE
     std::cout << "success: " << success << std::endl;
     std::cout << "result: " << res << std::endl;
     std::cout << "objective: " << objValue << std::endl;
     std::cout << "primal.size(): " << primal.size() << std::endl;
     std::cout << "dual.size(): " << dual.size() << std::endl;
#endif
     for(long i=0; i< primal.size(); i++){
       primal[i] = ( (double) simplex.flow( graph.arcFromId( i ) ) ) / capacityScaling;
       //primal[i] = ns.flow( graph.arcFromId( i ) );
     }
     for(long i=0; i<dual.size(); i++){
       dual[i] = simplex.potential( graph.nodeFromId(i) );
     }


   };



   virtual bool isOptimal(){
     return success;
   };


   virtual double getObjectiveValue(){
     return objValue;
   };


   virtual long getIterationCount(){
     return iCount;
   };


   virtual long getNumberOfRows(){
     return dual.size();
   };


   virtual long getNumberOfColumns(){
     return primal.size();
   };


   virtual void setupStandardBasis(){
     for(long i= 0; i< getNumberOfColumns(); i++){
       colStatus[i] = LPSolver::LOWER;
     }
     for(long i= 0; i< getNumberOfRows(); i++){
       rowStatus[i] =  LPSolver::BASIC;
     }
   };




   virtual void createLP(long nSource, long nTarget){
     deleteLP();
     ns=nSource;
     nt=nTarget;

   };




   virtual double getColumnPrimal(long col){
     return primal[col];
   };


   virtual void setRowBounds(long i, double m){
     mass[i] = m;
   };


   virtual double getRowBounds(long i){
     return mass[i];
   };


   virtual void setColumnBoundsLower(long col, double lb){
     colLB[col] = lb;
     colUB[col] = std::numeric_limits<double>::max();
   };
   virtual void setColumnBounds(long col, double lb, double ub){
     colLB[col] = lb;
     colUB[col] = ub;
   };

   virtual Status getColumnStatus(long col){
     return  colStatus[col];

   };

   virtual void setColumnStatus(long col, Status s){
     colStatus[col] = s;
   };


   virtual void setRowStatus(long row, Status s){
     rowStatus[row] = s;
   };


   virtual void setColumnObjective(long i, double cost){
     coeff[i] = cost;
   };

   virtual void setColumnCoefficients(long col, long s, long t){
     tInd[col] = t;
     sInd[col] = s;
   };

   virtual long getColumn(long col, long *ind, double *val){
     ind[0] = sInd[col];
     ind[1] = tInd[col];
     val[0] = 1;
     val[1] = -1;
     return 2;
   };

   virtual void addColumns(long n){
     sInd.resize( sInd.size() + n, -1  );
     tInd.resize( tInd.size() + n, -1 );
     coeff.resize( coeff.size() + n, 0 );
     primal.resize( primal.size() + n, 0 );
     colStatus.resize( colStatus.size() + n, LPSolver::LOWER );
     colLB.resize( colLB.size() + n, 0 );
     colUB.resize( colUB.size() + n, 1 );
   };



   virtual void addRows(long n){
     mass.resize( mass.size() + n , 0);
     rowStatus.resize( rowStatus.size() + n, LPSolver::BASIC );
     dual.resize( dual.size() + n, 0);
   };



   virtual double getRowDual(long row){
     return dual[row];
   };

   virtual Status getRowStatus(long row){
     return rowStatus[row];
   };



  private:
   
   virtual void deleteLP(){
     sInd.clear();
     tInd.clear();
     coeff.clear();
     mass.clear();
     dual.clear();
     primal.clear();
     rowStatus.clear();
     colStatus.clear();
     colLB.clear();
     colUB.clear();
     ns = 0;
   };



};


#endif
