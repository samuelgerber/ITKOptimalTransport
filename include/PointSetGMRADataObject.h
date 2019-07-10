#ifndef POINTSETGMRADATAOBJECT_H
#define POINTSETGMRADATAOBJECT_H


#include <Eigen/Dense>



template <typename TPointSetType>
class PointSetGMRADataObject : public GMRADataObject<typename TPointSetType::PixelType>{
  public:

    using PointSetType = TPointSetType;
    using PointType = typename TPointSetType::PointType;
    using PointSetPointer = typename TPointSetType::Pointer;
    using PixelType = typename TPointSetType::PixelType;

    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<PixelType, Eigen::Dynamic, 1> VectorXp;

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

