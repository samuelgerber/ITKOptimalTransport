#ifndef TRANSPORTLPSOLVER_H
#define TRANSPORTLPSOLVER_H


#include "MultiscaleTransport.h"
#include "LPSolver.h"


template <typename TPrecision>
class TransportLPSolver {

   public:

    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;
    typedef typename TransportPlan<TPrecision>::Path Path;
  
    typedef LPSolver::Status Status; 

    enum TransportType { BALANCED,
                         UNBALANCED_ADD,
                         UNBALANCED_SUBTRACT,
                         UNBALANCED_SOURCE,
                         UNBALANCED_FREE,
                         UNBALANCED_FREE_SOURCE
                       };




  private:
    
    int pathOffset;

    int sourceMassSupplyDistributionStart;
    int sourceMassSinkDistributionStart;
    int targetMassSupplyDistributionStart;
    int targetMassSinkDistributionStart;
    //int sourceMassSupplyPath;
    //int sourceMassSinkPath;
    //int targetMassSupplyPath;
    //int targetMassSinkPath;
    int sourceMassExchangePath;
    int targetMassExchangePath; 

    int sourceCirculationPath;


    //int sourceTerminalNode;
    int sourceMassSupplyNode;
    int sourceMassSinkNode;
    //int targetTerminalNode;
    int targetMassSupplyNode;
    int targetMassSinkNode;

    TransportType transportType;
    double massDeltaCost;
    double lambda;

    LPSolver *solver;
    bool lastScale;

  public:

   TransportLPSolver(LPSolver *lpSolver, TransportType type, double massCost = 0, double l=0) :
      solver(lpSolver),
      transportType(type),
      massDeltaCost( massCost ),
      lambda(l),
      lastScale(false)
   {

   };

   virtual ~TransportLPSolver(){};

   void setLastScale(bool last){
     lastScale=last;
   };

   void createLP( TransportPlan<TPrecision> *sol ){

#ifdef VERBOSE
      std::cout << "TransportLPSolver Setup" << std::endl;
      std::cout << " lambda: " << lambda << std::endl;
      std::cout << " massCost: " << massDeltaCost << std::endl;
      std::cout << " TransportType: " << transportType << std::endl;
#endif
      int ns = sol->source->getNodes().size();
      int nt = sol->target->getNodes().size();

      //Source and Target masses
      double massSource = 0;
      double massTarget = 0;
      for( TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
           sol->source->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        massSource += n->getMass();
      }
      for( TransportNodeVectorCIterator it = sol->target->getNodes().begin(); it !=
           sol->target->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        massTarget += n->getMass();
      }
      double massDelta = fabs(massSource - massTarget);
      double massTotal = std::max( massSource, massTarget);
      double massMinimal = std::min( massSource, massTarget);

#ifdef VERBOSE
      std::cout << " massTarget: " << massTarget << std::endl;
      std::cout << " massSource: " << massSource << std::endl;
      std::cout << " massDelta: " << massDelta << std::endl;
#endif

      //Setup LP
      solver->createLP( ns + 3, nt + 3 );
      
      //rows - constraints:
      int nConstraints = ns + nt + 4;
      solver->addRows( nConstraints );
      
      //Node IDs:
      // - Sources: First ns nodes
      // - Targets: Next nt nodes
      // - Source/Target terminal and unbalanced nodes:
      //sourceTerminalNode   = ns+nt;
      sourceMassSupplyNode = ns+nt+0;
      sourceMassSinkNode   = ns+nt+1;
      //targetTerminalNode   = ns+nt+3;
      targetMassSupplyNode = ns+nt+2;
      targetMassSinkNode   = ns+nt+3;



      //columns - paths between source and targets plus paths to terminal nodes
      //and 2 mass paths for each terminal node 
      solver->addColumns( sol->getNumberOfPaths() + 2*ns + 2*nt + 3 ); 
      //Path ids 
      sourceMassSupplyDistributionStart = 0;
      sourceMassSinkDistributionStart = ns;
      targetMassSupplyDistributionStart = 2*ns;
      targetMassSinkDistributionStart = 2*ns + nt;
      int distributionEnd = targetMassSupplyDistributionStart + 2*nt;
      //sourceMassSupplyPath    = targetMassSupplyDistributionStart + 2*nt;
      //sourceMassSinkPath      = sourceMassSupplyPath + 1;
      //targetMassSupplyPath    = sourceMassSupplyPath + 2;
      //targetMassSinkPath      = sourceMassSupplyPath + 3;
      sourceMassExchangePath  = distributionEnd;
      targetMassExchangePath  = distributionEnd + 1;
      sourceCirculationPath  = distributionEnd + 2;
      pathOffset              = distributionEnd + 3;
      
      //Mass suply and distribution upper an lower bounds
      //double sourceMassSupplyLower = 0;
      //double sourceMassSupplyUpper = 0;
      //double sourceMassSinkLower = 0;
      //double sourceMassSinkUpper = 0;
      double sourceMassExchangeLower = 0;
      double sourceMassExchangeUpper = 0;
      double sourceMassSupply = 0;
      double sourceMassSink = 0;

      double sourceCirculationUpper = 0;
      
      //double targetMassSupplyLower = 0;
      //double targetMassSupplyUpper = 0;
      //double targetMassSinkLower = 0;
      //double targetMassSinkUpper = 0;
      double targetMassExchangeLower = 0;
      double targetMassExchangeUpper = 0;
      double targetMassSupply = 0;
      double targetMassSink = 0;

      //if(lastScale){
      if( transportType == UNBALANCED_ADD ){
        if( massSource > massTarget ){
          //targetMassSinkUpper = massDelta;
          //targetMassSinkLower = massDelta;
          targetMassSink = -massDelta;
        }
        else{
          //sourceMassSupplyUpper = massDelta;
          //sourceMassSupplyLower = massDelta;
          sourceMassSupply = massDelta;
        }
      }
      else if( transportType == UNBALANCED_SUBTRACT  ){
        if( massSource > massTarget ){
          //sourceMassSinkUpper = massDelta;
          //sourceMassSinkLower = massDelta;
          sourceMassSink = -massDelta;
        }
        else{
          //targetMassSupplyUpper = massDelta;
          //targetMassSupplyLower = massDelta;
          targetMassSupply = massDelta;
        }
      }
      else if( transportType == UNBALANCED_SOURCE ){
        if( massSource > massTarget ){
          //sourceMassSinkUpper = massDelta;
          //sourceMassSinkLower = massDelta;
          sourceMassSink = -massDelta;
        }
        else{
          //sourceMassSupplyUpper = massDelta;
          //sourceMassSupplyLower = massDelta;
          sourceMassSupply = massDelta;
        }
      }
      else if( transportType == UNBALANCED_FREE ){
        //sourceMassSupplyLower = 0;
        //sourceMassSupplyUpper = massTotal;
        sourceMassSupply = massTotal;
        //sourceMassSinkLower = 0;
        //sourceMassSinkUpper = massSource;
        sourceMassSink = -massSource + 0. * massMinimal;

        //targetMassSupplyLower = 0;
        //targetMassSupplyUpper = massTarget;
        targetMassSupply = massTarget - 0. * massMinimal;
        //targetMassSinkLower = 0;
        //targetMassSinkUpper = massTotal;
        targetMassSink = -massTotal;
        
        sourceMassExchangeLower = 0;
        sourceMassExchangeUpper = massTotal;
        targetMassExchangeLower = 0;
        targetMassExchangeUpper = massTotal;
      }
      else if( transportType == UNBALANCED_FREE_SOURCE){
        if( massSource > massTarget ){
          //sourceMassSinkUpper = massDelta;
          //sourceMassSinkLower = massDelta;
          sourceMassSink = -massDelta;
          //sourceMassSupplyUpper = massDelta;
          //sourceMassSupplyLower = 0;

        }
        else{
          //sourceMassSinkUpper = massDelta;
          //sourceMassSinkLower = 0;
          //sourceMassSupplyUpper = massDelta;
          //sourceMassSupplyLower = massDelta;
          sourceMassSupply = massDelta;
        }
        sourceCirculationUpper = massSource;

      }
      //}


      //Source supply and terminal nodes and paths
      //solver->setRowBounds( sourceTerminalNode, 0);//massSource);
      solver->setRowBounds( sourceMassSupplyNode,  sourceMassSupply );
      solver->setRowBounds( sourceMassSinkNode, sourceMassSink      );

      //solver->setColumnBounds( sourceMassSupplyPath, 
      //    sourceMassSupplyLower, sourceMassSupplyUpper );
      //solver->setColumnObjective( sourceMassSupplyPath,  0);//massDeltaCost );
      //solver->setColumnCoefficients( sourceMassSupplyPath, 
      //    sourceMassSupplyNode, sourceTerminalNode);//, 1, -1 );

      //solver->setColumnBounds( sourceMassSinkPath, 
      //    sourceMassSinkLower, sourceMassSinkUpper );
      //solver->setColumnObjective( sourceMassSinkPath,  0);//massDeltaCost );
      //solver->setColumnCoefficients( sourceMassSinkPath, 
      //    sourceTerminalNode, sourceMassSinkNode);//, 1, -1 );



      //Target supply and terminal nodes and paths
      //solver->setRowBounds( targetTerminalNode, 0);//-massTarget);
      solver->setRowBounds( targetMassSupplyNode, targetMassSupply );
      solver->setRowBounds( targetMassSinkNode,   targetMassSink   );

      //solver->setColumnBounds( targetMassSupplyPath,
      //    targetMassSupplyLower, targetMassSupplyUpper );
      //solver->setColumnObjective( targetMassSupplyPath,  massDeltaCost );
      //solver->setColumnCoefficients( targetMassSupplyPath, 
      //    targetMassSupplyNode, targetTerminalNode);//, 1, -1 );

      //solver->setColumnBounds( targetMassSinkPath, 
      //    targetMassSinkLower, targetMassSinkUpper );
      //solver->setColumnObjective( targetMassSinkPath,  0);//massDeltaCost );
      //solver->setColumnCoefficients( targetMassSinkPath, 
      //    targetTerminalNode, targetMassSinkNode );//, 1, -1 );


      //Mass exchange paths
      solver->setColumnBounds( sourceMassExchangePath, 
          sourceMassExchangeLower, sourceMassExchangeUpper);
      solver->setColumnObjective( sourceMassExchangePath,  0 );
      solver->setColumnCoefficients( sourceMassExchangePath, 
          sourceMassSupplyNode, targetMassSinkNode);//, 1, -1 );


      solver->setColumnBounds( targetMassExchangePath, 
          targetMassExchangeLower, targetMassExchangeUpper);
      solver->setColumnObjective( targetMassExchangePath,  0 );
      solver->setColumnCoefficients( targetMassExchangePath, 
          targetMassSupplyNode, sourceMassSinkNode );//, 1, -1 );

      //Circulation path for redistribution mass in source
      solver->setColumnBounds( sourceCirculationPath, 0, sourceCirculationUpper);
      solver->setColumnObjective( sourceCirculationPath,  0);
      solver->setColumnCoefficients( sourceCirculationPath, 
          sourceMassSinkNode, sourceMassSupplyNode );//, 1, -1 );


      //Row constraints amount of mass in / out flow at each node
      double wfactor = 1 + 10 * massDelta / massTotal;
      for(TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
          sol->source->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        double w = n->getMass();
        double wdelta = lambda*w;
        double lower = 0;//std::max(0., w  - wdelta - sourceMassSinkUpper);
        double upper = w;// + wdelta + sourceMassSupplyUpper + targetMassSinkUpper;
        
        solver->setRowBounds(n->getID(),  w);

        int pid = sourceMassSupplyDistributionStart + n->getID();
        solver->setColumnBoundsLower( pid, lower);
        solver->setColumnBounds( pid, lower, 0.99999 * upper);
        solver->setColumnObjective( pid,  massDeltaCost );
        solver->setColumnCoefficients( pid, sourceMassSupplyNode, n->getID() ); //, 1, -1 );

        pid = sourceMassSinkDistributionStart + n->getID();
        solver->setColumnBoundsLower( pid, lower);
        solver->setColumnBounds( pid, lower, 0.99999*upper);
        solver->setColumnObjective( pid,  massDeltaCost );
        solver->setColumnCoefficients( pid, n->getID(), sourceMassSinkNode); //, 1, -1 );
      }

      int offset= sol->source->getNodes().size();

      for(TransportNodeVectorCIterator it = sol->target->getNodes().begin(); it !=
          sol->target->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        double w = n->getMass();
        double wdelta = lambda*w;
        double lower = 0;//std::max(0., w - wdelta - targetMassSupplyUpper);
        double upper = w;// + wdelta +  sourceMassSupplyUpper + targetMassSinkUpper;

        solver->setRowBounds( offset + n->getID(), -w);
        
        int pid = targetMassSinkDistributionStart + n->getID();
        solver->setColumnBoundsLower( pid, lower);//
        solver->setColumnBounds( pid, lower, 0.99999 * upper);
        solver->setColumnObjective( pid,  massDeltaCost );
        solver->setColumnCoefficients( pid, offset + n->getID(), targetMassSinkNode ); //, 1, -1  );

        pid = targetMassSupplyDistributionStart + n->getID();
        solver->setColumnBoundsLower( pid, lower);
        solver->setColumnBounds( pid, lower, 0.99999*upper);
        solver->setColumnObjective( pid,  massDeltaCost );
        solver->setColumnCoefficients( pid, targetMassSupplyNode, offset + n->getID() ); //, 1, -1  );
      }


      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        Path &path = sol->pathIteratorCurrent();
        TransportNode<TPrecision> *from = path.from;
        TransportNode<TPrecision> *to = path.to;
        this->setColumnBoundsLower(path.index, 0);
        this->setColumnObjective(path.index, path.cost );
        this->setColumnCoefficients(path.index, from->getID(), offset + to->getID());
      }



   }




   void storeLP(TransportPlan<TPrecision> *sol, TPrecision p){
     sol->cost = pow( this->getObjectiveValue(), 1.0/p );

     //Store solution
     double sumw = 0;
     int nNonZero =0;
     for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
         sol->pathIteratorNext() ){
       Path &path = sol->pathIteratorCurrent();
       path.w = this->getColumnPrimal(path.index);
       sumw += path.w;
       nNonZero += (path.w > 0);
     }

#ifdef VERBOSE
     std::cout << "cost: " << sol->cost << std::endl;
     std::cout << "nonzeros: "      << nNonZero << std::endl;
     std::cout << "solution sumw: " << sumw     << std::endl << std::endl;
#endif

     //Potential of from nodes
     for(TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
         sol->source->getNodes().end(); ++it ){
       TransportNode<TPrecision> *node = *it;
       TPrecision pi = this->getRowDual( node->getID() );
       node->setPotential(pi);
     }


     //Potential of to nodes
     int offset = sol->source->getNodes().size();
     for(TransportNodeVectorCIterator it = sol->target->getNodes().begin(); it !=
         sol->target->getNodes().end(); ++it){
       TransportNode<TPrecision> *node = *it;
       TPrecision pi = this->getRowDual( offset + node->getID() );
       node->setPotential(pi);
     }

     
/*
     double massAllocation = 0;     
     double massDeallocation = 0;     
     for(TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
          sol->source->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;

        int pid = sourceMassSupplyDistributionStart + n->getID();
        massAllocation += solver->getColumnPrimal(pid);

        pid = sourceMassSinkDistributionStart + n->getID();
        massDeallocation += solver->getColumnPrimal(pid);
     }

     std::cout << "Allocation: " << massAllocation << std::endl;
     std::cout << "Deallocation: " << massDeallocation << std::endl;
*/

   };


   void setPotentials( TransportNodeVector &sourceNodes, 
                              TransportNodeVector &targetNodes){

      for(TransportNodeVectorIterator sIt = sourceNodes.begin(); sIt !=
          sourceNodes.end(); ++sIt){
        TransportNode<TPrecision> *node = *sIt;
        TPrecision pi = this->getRowDual(node->getID());
        node->setPotential(pi);
        node->resetPi(pi, pi);
      }

      int offset = sourceNodes.size();
      for(TransportNodeVectorIterator tIt = targetNodes.begin(); tIt !=
          targetNodes.end(); ++tIt){
        TransportNode<TPrecision> *node = *tIt;
        TPrecision pi = this->getRowDual( offset + node->getID() );
        node->setPotential(pi);
        node->resetPi(pi, pi);
      }

    };




    void setupBasis(std::vector<Status> &colStatus, std::vector<Status> &rowStatus){

      //int nBasic = 0;
      for(int i=0; i<colStatus.size(); i++){
        solver->setColumnStatus(i, colStatus[i]);
        //nBasic += colStatus[i] == LPSolver::BASIC;
      }

      for(int i=0; i<rowStatus.size(); i++){
        solver->setRowStatus(i, rowStatus[i]);
        //nBasic += rowStatus[i] == LPSolver::BASIC;
      }


    };

   

    void setColumnBoundsLower(long col, double lb){
      solver->setColumnBoundsLower( pathOffset + col, lb);
    };

    void setColumnBounds(long col, double lb, double ub){
      solver->setColumnBounds( pathOffset + col, lb, ub);
    };

    void setColumnObjective(long col, double cost){
      solver->setColumnObjective(pathOffset + col, cost);
    };


    void setColumnCoefficients( long col, long s, long t){
      solver->setColumnCoefficients(pathOffset + col, s, t);
    
    };

    double getColumnPrimal(long col){
      return solver->getColumnPrimal( pathOffset + col );
    };


    Status getColumnStatus(long col){
      return solver->getColumnStatus(pathOffset + col);
    };

    void setColumnStatus(long col, Status s){
      solver->setColumnStatus(pathOffset+col, s);
    };

    long getColumn(long col, long *ind, double *val){
      return solver->getColumn(pathOffset + col, ind, val);
    }

    void addColumns(long n){
      solver->addColumns(n);
    }

    Status getRowStatus(long row){
      return solver->getRowStatus(row);
    }

    void setRowStatus(long row, Status s){
      solver->setRowStatus(row, s);
    }
 
    double getRowDual(long row){
      return solver->getRowDual(row);
    }

    void solveLP(){
      solver->solveLP();
#ifdef VERBOSE
      std::cout << "sourceMassExchnage: " << solver->getColumnPrimal( sourceMassExchangePath ) << std::endl;
      std::cout << "targetMassExchange: " << solver->getColumnPrimal( targetMassExchangePath ) << std::endl;
      std::cout << "Source Circulation: " << solver->getColumnPrimal( sourceCirculationPath ) << std::endl;
#endif


    }

    bool isOptimal(){
      return solver->isOptimal();
    };

    long getNumberOfRows(){
      return solver->getNumberOfRows();
    };

    long getNumberOfColumns(){
      return solver->getNumberOfColumns();
    };

    double getObjectiveValue(){
      return solver->getObjectiveValue();
    };

    long getIterationCount(){
      return solver->getIterationCount();
    };
  
    void setupStandardBasis(){
      solver->setupStandardBasis();
    };
   
};


#endif


