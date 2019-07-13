#ifndef TRANSPORTPLAN_H 
#define TRANSPORTPLAN_H 

#include "TransportNode.h"
#include "MultiscaleTransportLevel.h"

#include <map>
#include <vector>
#include <list>

#include <Eigen/Dense>

template <typename TPrecision>
class TransportPlan{



  public:


    struct Path{
      TransportNode<TPrecision> *from;
      TransportNode<TPrecision> *to;

      Path() {
        from = NULL;
        to = NULL;
        index= -1;
        cost = -1;
        w = 0;
      };


      Path(TransportNode<TPrecision> *f, TransportNode<TPrecision> *t)
        :from(f),to(t){
          index= -1;
          cost = -1;
          w = 0;
        };


      bool operator < (const Path &s1) const{
        if(from->getID() < s1.from->getID() ){
          return true;
        }
        else if(from->getID() > s1.from->getID() ){
          return false;
        }
        else if(to->getID() < s1.to->getID() ){
          return true;
        }
        return false;
      };


      bool operator > (const Path &s1) const {
        if(from->getID() > s1.from->getID() ){
          return true;
        }
        else if(from->getID() < s1.from->getID() ){
          return false;
        }
        else if(to->getID() > s1.to->getID() ){
          return true;
        }
        return false;
      };

      bool operator == (const Path &s1) const{
        return from->getID() == s1.from->getID() && to->getID() == s1.to->getID();
      };


      int index;
      TPrecision cost;
      TPrecision w;

    };




  private:

    std::vector< std::map< int, Path> > paths;
    typename std::vector< std::map<int, Path> >::iterator pathIteratorOuter;
    typename std::map< int, Path>::iterator pathIteratorInner;

    std::vector<int> toPathCounts;
    int pathCounter;



  public:


    TransportPlan(MultiscaleTransportLevel<TPrecision> *s,
        MultiscaleTransportLevel<TPrecision> *t) : source(s), target(t){

      paths.resize( source->getNodes().size() );
      toPathCounts.resize( target->getNodes().size(), 0 );
      timeSolve = 0;
      timeRefine = 0;
      timePropagate = 0;
      cost = std::numeric_limits<TPrecision>::max();
      optimizationStatus = -1;
      //nTotalPaths = -1;
      pathCounter = 0;

    };

/*
    TransportPlan(int nFrom, int nTo) : source(NULL), target(NULL){

      paths.resize( nFrom );
      toPathCounts.resize( nTo, 0);
      timeLP = 0;
      timeNeighborhood = 0;
      timePropagate = 0;
      cost = std::numeric_limits<TPrecision>::max();
      optimizationStatus = -1;
      nTotalPaths = -1;
      pathCounter = 0;
    };

*/

    virtual ~TransportPlan(){
    };


    MultiscaleTransportLevel<TPrecision> *source;
    MultiscaleTransportLevel<TPrecision> *target;


    double cost;
    int optimizationStatus;
    //int nTotalPaths;
    clock_t timeSolve;
    clock_t timePropagate;
    clock_t timeRefine;


    //For sinkhorn transport
    Eigen::MatrixXd leftScaling;


    bool hasPath(const Path &p){
       return hasPath( p.from->getID(), p.to->getID() );
    };


    bool hasPath(int from, int to){
      return paths[from].find(to) != paths[from].end();
    };


    Path &getPath(TransportNode<TPrecision> *from, TransportNode<TPrecision> *to){
       return getPath( from->getID(), to->getID() );
    };


    Path &getPath(int from, int to){
      return paths[from][to];
    };


    int getNumberOfToPaths(int from){
      return paths[from].size();
    };
    int getNumberOfFromPaths(int to){
      return toPathCounts[to];
    };

    int addPath(Path p){
      if( !hasPath(p) ){
        p.index = pathCounter;
        ++pathCounter;
        paths[p.from->getID()][p.to->getID()] = p;
        toPathCounts[p.to->getID()] += 1;
        return p.index;
      }
      else{
        Path &path = getPath(p.from, p.to);
        path.w = std::max(path.w, p.w);
        return path.index;
      }
    };


    int getPathIndex(const Path &p){
      int from = p.from->getID();
      int to = p.to->getID();

      typename std::map<int, Path>::const_iterator it = paths[from].find(to);

      if(it == paths[from].end()){
        return -1;
      }
      return it->second.index;
    };


    void pathIteratorBegin(){

      pathIteratorOuter = paths.begin();
      if(pathIteratorOuter != paths.end() ){

        pathIteratorInner = (*pathIteratorOuter).begin();

        while( pathIteratorInner == (*pathIteratorOuter).end() ){
          ++pathIteratorOuter;
          if( pathIteratorOuter != paths.end()){
            pathIteratorInner = (*pathIteratorOuter).begin();
          }
          else{
            break;
          }
        }
      }

    };


    bool pathIteratorIsAtEnd(){
      return pathIteratorOuter == paths.end();
    };


    Path &pathIteratorCurrent(){
      return pathIteratorInner->second;
    };




    void pathIteratorNext(bool erase=false){
      if(erase){
        (*pathIteratorOuter).erase(pathIteratorInner++);
        pathCounter--;
      }
      else{
        ++pathIteratorInner;
      }
      while( pathIteratorInner == (*pathIteratorOuter).end() ){
        ++pathIteratorOuter;
        if( pathIteratorOuter != paths.end()){
          pathIteratorInner = (*pathIteratorOuter).begin();
        }
        else{
          break;
        }
      }
    };



    int getNumberOfPaths(){
      return pathCounter;
    };




    //Compute a multicsale path cost
    std::vector<TPrecision> getMultiscaleTransportCost(int p){
      std::vector<TPrecision> costs( std::max(source->getScale(), target->getScale())+1, 0 );

      for( this->pathIteratorBegin(); !this->pathIteratorIsAtEnd();
           this->pathIteratorNext() ){
        Path &path = this->pathIteratorCurrent();

        if(path.w > 0 ){
          typename std::vector< TransportNode<TPrecision>* > fNodes;
          TransportNode<TPrecision> *from = path.from;
          while(from != NULL){
            fNodes.push_back(from);
            from = from->getParent();
          }
          typename std::vector< TransportNode<TPrecision>* > tNodes;
          TransportNode<TPrecision> *to = path.to;
          while(to != NULL){
            tNodes.push_back(to);
            to = to->getParent();
          }


          typename std::vector< TransportNode<TPrecision>* >::reverse_iterator fIt = fNodes.rbegin();
          typename std::vector< TransportNode<TPrecision>* >::reverse_iterator tIt = tNodes.rbegin();
          from = *fIt;
          to = *tIt;

          int index = 0;
          while(tIt != tNodes.rend() || fIt != fNodes.rend()){
            costs[index] = costs[index] + from->getTransportCost(to, p) * path.w ;
            ++index;

            if(tIt != tNodes.rend()){
              ++tIt;
            }
            if(tIt != tNodes.rend() ){
              to = *tIt;
            }

            if(fIt != fNodes.rend()){
              ++fIt;
            }
            if(fIt != fNodes.rend() ){
              from = *fIt;
            }

          }

        }

      }

      return costs;
    };


    TransportPlan<TPrecision> *createCopy(){
      TransportPlan<TPrecision> *res = new
        TransportPlan<TPrecision>(source, target);

      res->timePropagate = this->timePropagate;
      res->timeSolve = this->timeSolve;
      res->timeRefine = this->timeRefine;
      res->paths = this->paths;
      res->pathCounter = this->pathCounter;
      res->toPathCounts = this->toPathCounts;
      return res;
    };

    void copyFrom(TransportPlan<TPrecision> *res){

      this->timePropagate = res->timePropagate;
      this->timeSolve = res->timeSolve;
      this->timeRefine = res->timeRefine;
      this->paths = res->paths;
      this->pathCounter = res->pathCounter;
      this->toPathCounts = res->toPathCounts;
      this->nTotalPaths = res->nTotalPaths;
      return res;
    };

};
















template <typename TPrecision>
class TransportPlanSolutions{

  private:

    typedef typename  std::list< TransportPlan<TPrecision> *>::iterator TPIterator;
    typedef typename TransportPlan<TPrecision>::Path Path;


    std::list< TransportPlan<TPrecision> *> alternatives;
    TransportPlan<TPrecision> *sol;



  public:




    TransportPlanSolutions(MultiscaleTransportLevel<TPrecision> *s,
        MultiscaleTransportLevel<TPrecision> *t) {
      sol = new TransportPlan<TPrecision>(s, t);
    };


    virtual ~TransportPlanSolutions(){

      for(TPIterator it = alternatives.begin(); it !=alternatives.end(); ++it){
        delete *it;
      }
      alternatives.clear();

    };


    void setPrimarySolution(TransportPlan<TPrecision> *primary){
      sol = primary;
    };

    TransportPlan<TPrecision> *getPrimarySolution(){
      return sol;
    };


    std::list< TransportPlan<TPrecision> * > &getAlternativeSolutions(){
      return alternatives;
    };


    void addAlternativeSolution( TransportPlan<TPrecision> *a){
      alternatives.push_back(a);
    };



    TransportPlan<TPrecision> *getCombinedPaths(){

      TransportPlan<TPrecision> *res = sol->createCopy();

      for(TPIterator it = alternatives.begin(); it !=alternatives.end(); ++it){
        TransportPlan<TPrecision> *a= *it;
        for( a->pathIteratorBegin(); !a->pathIteratorIsAtEnd();
                 a->pathIteratorNext() ){
          Path &path = a->pathIteratorCurrent();
          res->addPath(path);
        }
      }

      return res;
    };

};





#endif
