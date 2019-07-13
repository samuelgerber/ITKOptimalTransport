#ifndef POINTSETGMRADATAOBJECT_H
#define POINTSETGMRADATAOBJECT_H


#include <Eigen/Dense>
#include "GMRADataObject.h"


template <typename TPointSetType>
class PointSetGMRADataObject : public GMRADataObject<double>{
  public:

    using PointSetType = TPointSetType;
    using PointType = typename TPointSetType::PointType;
    using PointSetPointer = typename TPointSetType::Pointer;
    using CoordRepType = typename TPointSetType::CoordRepType;

    typedef typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXp;

    PointSetGMRADataObject(const PointSetType *ps) : pointSet(ps){};

    virtual VectorXp getPoint(int i){
       VectorXp p( this->dimension() );
       PointType tmp = pointSet->GetPoint(i);
       for(int i=0; i<this->dimension(); i++){
         p[i] = tmp[i];
       }
       return p;
    };

    virtual int numberOfPoints(){
      return pointSet->GetNumberOfPoints();
    };

    virtual int dimension(){
      return PointSetType::PointDimension;
    }

  private:
    const PointSetType *pointSet;

};

#endif

