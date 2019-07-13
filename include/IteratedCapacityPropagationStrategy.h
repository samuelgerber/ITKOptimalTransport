#ifndef ITERATEDCAPACITYPROPAGATIONSTRATEGY_H
#define ITERATEDCAPACITYPROPAGATIONSTRATEGY_H

#include "NeighborhoodPropagationStrategy.h"
#include "Random.h"

#include <ctime>


template <typename TPrecision>
class IteratedCapacityPropagationStrategy : public NeighborhoodPropagationStrategy<TPrecision>{



  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path;


  private:
    int nIterations;

  public:

    IteratedCapacityPropagationStrategy(int nIter, TPrecision eFactor=0) :
      NeighborhoodPropagationStrategy<TPrecision>(eFactor), nIterations(nIter){
    };

    virtual ~IteratedCapacityPropagationStrategy(){};




  protected:

    virtual void computeAlternateSolutions( TransportLPSolver<TPrecision> *solver,
                        TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale ){

      TransportPlan<TPrecision> *sol = pSol->getPrimarySolution();

      TransportPlan<TPrecision> *res = sol->createCopy();
      solver->createLP(res);

      for(int i=0; i<nIterations; i++){
        clock_t t1 = clock();
        for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
          res->pathIteratorNext() ){
          Path &path = res->pathIteratorCurrent();

          if(path.w > 0){
            static Random<TPrecision> random;
            int nFrom = res->getNumberOfToPaths( path.from->getID() );
            int nTo = res->getNumberOfFromPaths( path.to->getID() );
            int div = std::min(nTo, nFrom);
            TPrecision ub =  path.w * ( 0.9 + 0.05*random.Uniform() );
            //TPrecision ub =  path.w * ( 0.9 + 0.1 / (div-0.0001) );
            solver->setColumnBounds(path.index, 0, ub);
            //path.cost *= 10;
            //solver->setCoefficent(path.index, path.cost);
          }
        }
#ifdef VERBOSE
        std::cout << "Solving iterated capacity propagation problem" << std::endl;
#endif
        clock_t t2 = clock();
        solver->solveLP();
        clock_t t3 = clock();



        //Fix if unfeasible LP
        if( !solver->isOptimal() ){
#ifdef VERBOSE
          std::cout << "Suboptimal, problem infeasible. Removing some upper bounds" << std::endl;
#endif
          /*
          TransportNodeVector &snodes = res->source->getNodes();
          std::vector<double> solSNodeWeights(snodes.size(), 0);

          TransportNodeVector &tnodes = res->target->getNodes();
          std::vector<double> solTNodeWeights(tnodes.size(), 0);

          for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
               res->pathIteratorNext() ){
            Path &path = res->pathIteratorCurrent();

            if(path.w > 0){
              int fromID =  path.from->getID();
              solSNodeWeights[fromID] += path.w;
              int toID =  path.to->getID();
              solTNodeWeights[toID] += path.w;
            }
          }

          std::vector<double> reqSNodeWeights(snodes.size(), 0);
          for( TransportNodeVectorCIterator it = snodes.begin(); it != snodes.end();
                ++it ){
            reqSNodeWeights[ (*it)->getID() ] = (*it)->getMass();
          }

          std::vector<double> reqTNodeWeights(tnodes.size(), 0);
          for( TransportNodeVectorCIterator it = tnodes.begin(); it != tnodes.end();
                ++it ){
            reqTNodeWeights[ (*it)->getID() ] = (*it)->getMass();
          }


          int nBounds = 0;
          for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
               res->pathIteratorNext() ){
            Path &path = res->pathIteratorCurrent();

            if(path.w > 0){
              int fromID =  path.from->getID();
              int toID =  path.to->getID();
              if( solTNodeWeights[toID] < reqTNodeWeights[toID] ||
                  solSNodeWeights[fromID] < reqSNodeWeights[fromID] ){
	               solver->setColBounds(path.index);
                 nBounds++;
              }
            }
          }
          std::cout << "Removed " << nBounds << " upper bounds" << std::endl;
          std::cout << "Resolving" << std::endl;
          solver->solveLP();




          //reset bounds and recompute solution
          */
           for( res->pathIteratorBegin(); !res->pathIteratorIsAtEnd();
		            res->pathIteratorNext() ){
	          Path &path = res->pathIteratorCurrent();
	          solver->setColumnBoundsLower(path.index, 0);
          }
          solver->solveLP();

          //break;
        }

        solver->storeLP(res, p);
        pSol->addAlternativeSolution(res);
        res = res->createCopy();

      }

      delete res;

      //pSol->addAlternativeSolution(res);


    };




};


#endif


