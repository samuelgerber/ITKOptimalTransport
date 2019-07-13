#ifndef MULTISCALETRANSPORTLEVEL_H
#define MULTISCALETRANSPORTLEVEL_H

#include "TransportNode.h" 


template <typename TPrecision>
class MultiscaleTransportLevel{

  public:

    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;


  private:

    int scale;
    TransportNodeVector nodes;

    MultiscaleTransportLevel<TPrecision> *parent;

  public:

    MultiscaleTransportLevel(int s, MultiscaleTransportLevel<TPrecision> *parentLevel) : scale(s), parent(parentLevel){
    };


    virtual ~MultiscaleTransportLevel(){
      for(TransportNodeVectorIterator it = nodes.begin(); it != nodes.end(); ++it){
        delete *it;
      }
      nodes.clear();
    };


    void addNode(TransportNode<TPrecision> *node){
      nodes.push_back(node);
    };


    TransportNodeVector &getNodes(){
      return nodes;
    };


    virtual TransportNodeVector getNeighborhood(TransportNode<TPrecision> *node,
        TPrecision eps) const = 0;

    TPrecision getMaximalRadius(){
      TPrecision radius = 0;
      TransportNodeVector &nodes = getNodes();
      for(TransportNodeVectorCIterator it = nodes.begin(); it != nodes.end();
          ++it){
        if(radius <  (*it)->getRadius() ){
          radius = (*it)->getRadius();
        }
      }
      return radius;

    };



    int getScale(){
      return scale;
    };



    MultiscaleTransportLevel<TPrecision> *getRootLevel(){
      MultiscaleTransportLevel<TPrecision> *current = this;
      MultiscaleTransportLevel<TPrecision> *parent = current->getParent();
      while(parent!=NULL){
        current=parent;
        parent = current->getParent();
      }

      return current;
    };



    MultiscaleTransportLevel<TPrecision> *getParent(){
      return parent;
    };

};



#endif
