#ifndef MULTISCALETRANSPORT_H
#define MULTISCALETRANSPORT_H

#include "MultiscaleTransportLevel.h"
#include "TransportPlan.h"



template <typename TPrecision>
class MultiscaleTransport{

  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;


    typedef typename TransportPlan<TPrecision>::Path Path;


    virtual ~MultiscaleTransport(){};


    std::vector< TransportPlan<TPrecision> * > solve(std::vector<
        MultiscaleTransportLevel<TPrecision>* > &aLevels,
        std::vector<MultiscaleTransportLevel<TPrecision>* > &bLevels,
        double p = 1, int nScales1 = -1, int nScales2 = -1, bool matchStartLevel =
        false, bool scaleMass = true){

      //compute distances at each tree level
      typename std::vector< TransportPlan<TPrecision> *> solutions;

      if(nScales1 < 0 || nScales1 >= (int) aLevels.size()){
        nScales1 = aLevels.size()-1;
      }

      if(nScales2 < 0 || nScales2 >= (int) bLevels.size()){
        nScales2 = bLevels.size()-1;
      }

      normalize(aLevels, scaleMass);
      normalize(bLevels, scaleMass);

      typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator itA =
        aLevels.begin();
      rootSource = *itA;

      typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator itB =
        bLevels.begin();
      rootTarget = *itB;




      int scaleStart1 = aLevels.size()-1-nScales1;
      int scaleStart2 = bLevels.size()-1-nScales2;
      std::advance( itA, scaleStart1 );
      std::advance( itB, scaleStart2 );

      if(matchStartLevel){

        TPrecision rA = getMeanRadius( *itA, p );
        TPrecision rB = getMeanRadius( *itB, p );
        TPrecision delta = rA-rB;
        if(delta > 0){
          ++itA;
          while(itA != aLevels.end() ){
            rA = getMeanRadius(*itA, p);
            TPrecision tmp = rA-rB;
            if( fabs(tmp) > fabs(delta) ){
              --itA;
              break;
            }
            else{
              delta = tmp;
            }
            ++itA;
          }

          if(itA == aLevels.end() ){
            --itA;
          }

        }
        else{
          ++itB;
          while(itB != bLevels.end() ){
            rB = getMeanRadius(*itB, p);
            TPrecision tmp = rA-rB;
            if( fabs(tmp) > fabs(delta) ){
              --itB;
              break;
            }
            else{
              delta = tmp;
            }
            ++itB;
          }

          if(itB == bLevels.end() ){
            --itB;
          }
        }


      }


      MultiscaleTransportLevel<TPrecision> *A = *itA;
      MultiscaleTransportLevel<TPrecision> *B = *itB;

      TransportPlanSolutions<TPrecision> *prevSol = NULL;
      int currentScale = 0;
      while( itA != aLevels.end() || itB != bLevels.end() ){

#ifdef VERBOSE
	std::cout << std::endl << std::endl << "Scale: " << currentScale;
        std::cout << std::endl << std::endl;
#endif

        ++currentScale;

	//Check if this is the last scale
        bool lastScale = false;
        if(itA == aLevels.end() ){
          typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itB;
          it++;
          lastScale = it == bLevels.end();
        }
        else if(itB == bLevels.end() ){
          typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itA;
          it++;
          lastScale = it == aLevels.end();
        }
        else{
         typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itA;
          it++;
          lastScale = it == aLevels.end();
          it = itB;
          it++;
          lastScale = lastScale && it == bLevels.end();
        }


        TransportPlanSolutions<TPrecision> *sol = solveLP(A, B, prevSol, p, lastScale);
        solutions.push_back( sol->getPrimarySolution() );
        if( prevSol != NULL){
          delete prevSol;
        }
        prevSol = sol;


        if(itA != aLevels.end()){
         ++itA;
        }
        if( itA != aLevels.end()){
          A = *itA;
        }

        if( itB != bLevels.end()){
          ++itB;
        }
        if( itB != bLevels.end()){
          B = *itB;
        }

      }



      return solutions;
    };






  protected:


    virtual TransportPlanSolutions<TPrecision>
      *solveLP(MultiscaleTransportLevel<TPrecision> *source,
          MultiscaleTransportLevel<TPrecision> *target, TransportPlanSolutions<TPrecision>
          *prevSol, double p, bool lastScale) = 0;


    MultiscaleTransportLevel<TPrecision> *getRootSource(){
      return rootSource;
    };

    MultiscaleTransportLevel<TPrecision> *getRootTarget(){
      return rootTarget;
    };



  private:

    MultiscaleTransportLevel<TPrecision> *rootSource;
    MultiscaleTransportLevel<TPrecision> *rootTarget;


    void normalize(std::vector< MultiscaleTransportLevel<TPrecision> * > &levels, bool scaleMass){

      //Normalize masses to one at each level and match up parent child mass
      //relations
      for(typename std::vector< MultiscaleTransportLevel<TPrecision> *
          >::reverse_iterator it =
        levels.rbegin(); it != levels.rend(); ++it){
        TransportNodeVector &n = (*it)->getNodes();
        for(TransportNodeVectorCIterator nIt = n.begin(); nIt != n.end(); nIt++){
          const TransportNodeVector &kids = (*nIt)->getChildren();
          TPrecision sum = 0;
          //if(kids.size() == 0 ){
          //  sum = (*nIt)->getMass();
          //}
          for(TransportNodeVectorCIterator kIt = kids.begin(); kIt !=kids.end();
              ++kIt){
            sum += (*kIt)->getMass();
          }
          //std::cout << (*nIt)->getMass() - sum << std::endl;
          (*nIt)->setMass(sum);
        }
      }

      if(scaleMass){
      TPrecision total = 0;
      TransportNodeVector &nodes = levels[0]->getNodes();
      for(TransportNodeVectorIterator nIt = nodes.begin(); nIt != nodes.end(); nIt++){
        total = (*nIt)->getMass();
      }
      int end=levels.size();
      if( !scaleMass ){
         end = end-1;
      }
      for(int i=0; i<end; i++){
        TransportNodeVector &n = levels[i]->getNodes();
        for(TransportNodeVectorIterator nIt = n.begin(); nIt != n.end(); nIt++){
          (*nIt)->setMass( (*nIt)->getMass() / total );
        }
      }
      }

    };



    TPrecision getMeanRadius(MultiscaleTransportLevel<TPrecision> *level, double p){
        TransportNodeVector &n = level->getNodes();
        TPrecision radius = 0;
        for(TransportNodeVectorCIterator nIt = n.begin(); nIt != n.end(); nIt++){
          radius += (*nIt)->getLocalNodeRadius();
        }
        return radius / n.size();
    };


};


#endif
