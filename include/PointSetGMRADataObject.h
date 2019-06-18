#ifndef POINTSETGMRADATAOBJECT_H
#define POINTSGMRADATAOBJECT_H


#include <Eigen/Dense>



template <typename TPointSetType>
class PointSetGMRADataObject : public GMRADataObject<TPointSetType::PixelType>{
  public:

    using PointSetPointer = typename TPointSetType::Pointer;
    using PointSet = typename TPointSetType;
    using PixelType = typename TPointSetType::PixelType;

    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, 1> VectorXp;

    virtual PointSetGMRADataObject(PointSetPointer &ps) : pointSet(ps){};

    virtual VectorXp getPoint(int i){
       VectorXp p( this->dimension() )
       pointSet.GetPointData(i, p.data()
       return p;
    };

    virtual int numberOfPoints(){
      return pointSet.GetNumberOfPoints();
    };

    virtual int dimension(){
      return pointSet::PointDimension;
    }

  private:
    PointSetPointer &pointSet;

};

#endif

