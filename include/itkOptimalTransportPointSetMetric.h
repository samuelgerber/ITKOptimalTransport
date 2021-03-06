/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkOptimalTransportPointSetMetric_h
#define itkOptimalTransportPointSetMetric_h

#include "itkPointSetToPointSetMetric2v4.h"
#include "itkTransportCoupling.h"

namespace itk
{
/** \class OptimalTransportPointSetMetric
 * \brief Computes the Euclidan distance metric between two point sets.
 *
 *  Given two point sets the Euclidean distance metric (i.e. ICP) is
 *  defined to be the aggregate of all shortest distances between all
 *  possible pairings of points between the two sets.
 *
 *  We only have to handle the individual point case as the parent
 *  class handles the aggregation.
 *
 *  Reference:
 *    PJ Besl and ND McKay, "A Method for Registration of 3-D Shapes",
 *    IEEE PAMI, Vol 14, No. 2, February 1992
 *
 * \ingroup ITKMetricsv4
 */
template<typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet,
  class TInternalComputationValueType = double>
class ITK_TEMPLATE_EXPORT OptimalTransportPointSetMetric:
  public PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(OptimalTransportPointSetMetric);

  /** Standard class type aliases. */
  using Self = OptimalTransportPointSetMetric;
  using Superclass = PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet,
    TInternalComputationValueType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimalTransportPointSetMetric, PointSetToPointSetMetric2v4 );

  /** Types transferred from the base class */
  using MeasureType = typename Superclass::MeasureType;
  using DerivativeType = typename Superclass::DerivativeType;
  using LocalDerivativeType = typename Superclass::LocalDerivativeType;
  using PointType = typename Superclass::PointType;
  using PixelType = typename Superclass::PixelType;
  using PointIdentifier = typename Superclass::PointIdentifier;
  using TransportCouplingType = TransportCoupling< PointIdentifier, PointIdentifier, double >;
  using TransportEntry = typename TransportCouplingType::TransportEntry;
  using TransportMap = typename TransportCouplingType::TransportMap;
  /**
   * Calculates the local metric value for a single point.
   */
  MeasureType GetLocalNeighborhoodValue( const PointIdentifier , const PixelType & pixel = 0 ) const override;

  /**
   * Calculates the local value and derivative for a single point.
   */
  void GetLocalNeighborhoodValueAndDerivative( const PointIdentifier ,
    MeasureType &, LocalDerivativeType &, const PixelType & pixel = 0 ) const override;

  itkSetObjectMacro(Coupling, TransportCouplingType);
  itkGetObjectMacro(Coupling, TransportCouplingType);

protected:
  OptimalTransportPointSetMetric() = default;
  ~OptimalTransportPointSetMetric() override = default;

  /** PrintSelf function */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  bool RequiresMovingPointsLocator() const override
    {
    return false;
    };

  bool RequiresFixedPointsLocator() const override
    {
    return false;
    }
private:
  typename TransportCouplingType::Pointer m_Coupling;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOptimalTransportPointSetMetric.hxx"
#endif

#endif
