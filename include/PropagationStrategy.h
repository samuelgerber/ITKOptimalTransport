#ifndef PROPAGATIONSTRATEGY_H
#define PROPAGATIONSTRATEGY_H

#include "MultiscaleTransportLevel.h"
#include "TransportPlan.h"
#include "TransportLPSolver.h"

#include <ctime>


template <typename TPrecision>
class PropagationStrategy {


  public:


    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename TransportPlan<TPrecision>::Path Path; 



  public:


    virtual ~PropagationStrategy(){};

    virtual TransportPlanSolutions<TPrecision> *propagate(TransportLPSolver<TPrecision> *solver,
        MultiscaleTransportLevel<TPrecision> *source,
        MultiscaleTransportLevel<TPrecision> *target,
        TransportPlanSolutions<TPrecision> *pSol, double p, bool lastScale) = 0;

};


#endif


