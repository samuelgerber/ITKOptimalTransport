#ifndef POINTSETGMRADATAOBJECT_H
#define POINTSETGMRADATAOBJECT_H


#include <Eigen/Dense>



template <typename TPointSetType>
class PointSetGMRADataObject : public GMRADataObject<typename TPointSetType::PixelType>{
  public:

    using PointSetType = TPointSetType;
    using PointSetPointer = typename TPointSetType::Pointer;
    using PixelType = typename TPointSetType::PixelType;

    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, 1> VectorXp;

    PointSetGMRADataObject(PointSetPointer ps) : pointSet(ps){};

    virtual VectorXp getPoint(int i){
       VectorXp p( this->dimension() );
       pointSet->GetPointData(i, p.data() );
       return p;
    };

    virtual int numberOfPoints(){
      return pointSet->GetNumberOfPoints();
    };

    virtual int dimension(){
      return PointSetType::PointDimension;
    }

  private:
    PointSetPointer pointSet;

};

#endif

