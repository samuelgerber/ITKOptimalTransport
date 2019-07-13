#ifndef TRANSPORTNODE_H 
#define TRANSPORTNODE_H 

#include <vector>


template <typename TPrecision>
class TransportNode{

  public:
    typedef typename std::vector< TransportNode<TPrecision>*  >  TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;

  private:

    TPrecision piMax;
    TPrecision piMin;
    TPrecision potential;

    TPrecision mass;

    TransportNodeVector kids;
    TransportNode<TPrecision> *parent;

    TPrecision radius;

    int id;

    int scale;


  public:

    TransportNode(int nodeID, int sca):id(nodeID), scale(sca){
      piMin = std::numeric_limits<TPrecision>::max();
      piMax = -piMin;
      mass = -1;
      parent = NULL;
      radius = -1;
    };

    virtual ~TransportNode(){
    };


    TPrecision getMass() const{
      return mass;
    };


    void setMass(TPrecision m){
      //std::cout << m << std::endl;
      mass = m;
    };


    void addChild(TransportNode *node){
      kids.push_back(node);
    };


    TransportNodeVector getChildren() const{
      return kids;
    };


    void setParent(TransportNode<TPrecision> *p){
      parent = p;
    };


    TransportNode<TPrecision> *getParent() const{
      return parent;
    };


    virtual TPrecision getNodeRadius() const = 0;
    virtual TPrecision getLocalNodeRadius() const = 0;
    virtual std::vector<int> &getPoints() = 0;

    virtual TPrecision getTransportCostRadius(double p){

      //NOTE: will *not* update if children are added after calling getRadius
      if(radius < 0){
        radius = 0;
        for(TransportNodeVectorIterator it = kids.begin(); it != kids.end(); ++it){
          TPrecision d = getTransportCost(*it, p);
          if(d > radius){
            radius = d;
          }
        }
      }

      return radius;
    };


    TPrecision getAverageTransportCost(const TransportNode<TPrecision> *other, double p){
      TPrecision sum = 0;
      for(int i=0; i<kids.size(); i++){
        for(int j=0; j<other->kids.size(); j++){
          sum += kids[i]->getTransportCost(other->kids[j], p);
        }
      }

      return sum / ( kids.size() * other->kids.size() );

    };

    virtual TPrecision getTransportCost(const TransportNode<TPrecision> *other, double p) const= 0;
    //virtual TPrecision getDeltaTransportCost(TransportNode<TPrecision> *other, double p) const= 0;

    /*
    virtual bool operator == (const TransportNode<TPrecision> &other) const = 0;
    virtual bool operator <  (const TransportNode<TPrecision> &other) const = 0;
    virtual bool operator >  (const TransportNode<TPrecision> &other) const = 0;
    */
    //virtual size_t hashKey() const = 0;



    void resetPi(double piMi = std::numeric_limits<TPrecision>::max(), double
        piMa = - std::numeric_limits<TPrecision>::max() ){
      piMin = piMi;
      piMax = piMa;
    };

    void setPiMin(TPrecision pi){
      if(piMin > pi){
        piMin = pi;
      }
    };

    void setPiMax(TPrecision pi){
      if(piMax < pi){
        piMax = pi;
      }
    };

    TPrecision getPiMin(){
      return piMin;
    };

    TPrecision getPiMax(){
      return piMax;
    };

    void setPotential(TPrecision p){
      potential = p;
    };

    TPrecision getPotential() const{
      return potential;
    };

    int getID() const{
      return id;
    };

    int getScale() const{
      return scale;
    };


};





template <typename TPrecision>
class TransportNodeEqual{
  public:
    bool operator()(const TransportNode<TPrecision> *s1, const
        TransportNode<TPrecision> *s2) const{
      return s1->operator == (*s2);
    };
};




template <typename TPrecision>
class TransportNodeLess{
  public:
    bool operator()(const TransportNode<TPrecision> *s1, const
        TransportNode<TPrecision> *s2) const{
      return s1->operator < (*s2);
    };
};




template <typename TPrecision>
class TransportNodeHash{
  public:
    bool operator()(const TransportNode<TPrecision> *s1) const{
      return s1->hashKey();
    };
};














#endif
