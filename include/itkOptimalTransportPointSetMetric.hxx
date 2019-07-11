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
#ifndef itkOptimalTransportPointSetMetric_hxx
#define itkOptimalTransportPointSetMetric_hxx

#include "itkOptimalTransportPointSetMetric.h"

namespace itk
{

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename OptimalTransportPointSetMetric<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::MeasureType
OptimalTransportPointSetMetric<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValue( const PointIdentifier fixedIndex, const PixelType & itkNotUsed( pixel ) ) const
{
  PointType fixedPoint = this->m_FixedTransformedPointSet->GetPoint(fixedIndex);

  TransportEntry &targets = m_Coupling->GetMap()[fixedIndex];
  MeasureType distance = 0;
  for(typename TransportEntry::iterator it = targets.begin(); it != targets.end(); ++it)
    {
    PointType movingPoint = this->m_MovingTransformedPointSet->GetPoint( it->first );
    distance += it->second * fixedPoint.EuclideanDistanceTo( movingPoint );
    }
  return distance;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
OptimalTransportPointSetMetric<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValueAndDerivative( const PointIdentifier fixedIndex,
  MeasureType &measure, LocalDerivativeType & localDerivative, const PixelType & itkNotUsed( pixel ) ) const
{

  PointType fixedPoint = this->m_FixedTransformedPointSet->GetPoint(fixedIndex);

  TransportEntry &targets = m_Coupling->GetMap()[fixedIndex];
  MeasureType distance = 0;
  PointType derivative;
  derivative.Fill(0);
  for(typename TransportEntry::iterator it = targets.begin(); it != targets.end(); ++it)
    {
    PointType movingPoint = this->m_MovingTransformedPointSet->GetPoint( it->first );
    distance += it->second * fixedPoint.EuclideanDistanceTo( movingPoint );
    derivative += (movingPoint - fixedPoint) * it->second;
    }
  localDerivative = derivative;
  measure = distance;
}

/** PrintSelf method */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
OptimalTransportPointSetMetric<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
