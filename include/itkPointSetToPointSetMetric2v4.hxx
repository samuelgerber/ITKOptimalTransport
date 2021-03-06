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
#ifndef itkPointSetToPointSetMetric2v4_hxx
#define itkPointSetToPointSetMetric2v4_hxx

#include "itkPointSetToPointSetMetric2v4.h"
#include "itkIdentityTransform.h"
#include "itkCompensatedSummation.h"

namespace itk
{

        /** Constructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PointSetToPointSetMetric2v4()
{
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::MeasureType
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValue() const
{
  this->InitializeForIteration();


  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->m_VirtualTransformedPointSet->GetNumberOfPoints() != this->m_FixedTransformedPointSet->GetNumberOfPoints() )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }
  /*
   * Split pointset in nWorkUnit ranges and sum individually
   * This splitting is required in order to avoid having the threads
   * repeatedly write to same location causing false sharing
   */
  //Use STL container to make sure no unesecarry checks are performed
  using FixedTransformedVectorContainer = typename FixedPointsContainer::STLContainerType;
  using VirtualPointsContainer = typename VirtualPointSetType::PointsContainer;
  using VirtualVectorContainer =  typename VirtualPointsContainer::STLContainerType;
  const VirtualVectorContainer &virtualTransformedPointSet =
    this->m_VirtualTransformedPointSet->GetPoints()->CastToSTLConstContainer();
  const FixedTransformedVectorContainer &fixedTransformedPointSet =
    this->m_FixedTransformedPointSet->GetPoints()->CastToSTLConstContainer();

  PointIdentifierRanges ranges = this->CreateRanges();
  std::vector< CompensatedSummation< MeasureType > > threadValues( ranges.size() );
  std::function< void(unsigned int) > sumNeighborhoodValues =
      [ this, &threadValues, &ranges, &virtualTransformedPointSet, &fixedTransformedPointSet]
      (unsigned int rangeIndex)
    {
    CompensatedSummation< MeasureType > threadValue = 0;
    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );


    for( PointIdentifier index=ranges[rangeIndex].first; index < ranges[rangeIndex].second; index++)
      {
      if( this->IsInsideVirtualDomain( virtualTransformedPointSet[index] ) )
        {
        if( this->m_UsePointSetData )
          {
          bool doesPointDataExist = this->m_FixedPointSet->GetPointData( index, &pixel );
          if( ! doesPointDataExist )
            {
            itkExceptionMacro( "The corresponding data for point (pointId = " << index << ") does not exist." );
            }
          }
        threadValue += this->GetLocalNeighborhoodValue( index, pixel );
        }
      }
      threadValues[rangeIndex] = threadValue;
    };

  //Sum per thread
  MultiThreaderBase::New()->ParallelizeArray( (PointIdentifier) 0,
                          (PointIdentifier) ranges.size(),
                           sumNeighborhoodValues, nullptr );
  //Join sums
  CompensatedSummation<MeasureType> value = 0;
  for(unsigned int i=0; i < threadValues.size(); i++)
    {
    value += threadValues[i];
    }

  DerivativeType derivative;
  MeasureType valueSum = value.GetSum();
  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    valueSum /= static_cast<MeasureType>( this->m_NumberOfValidPoints );
    }
  this->m_Value = valueSum;

  return valueSum;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetDerivative( DerivativeType & derivative ) const
{
  MeasureType value = NumericTraits<MeasureType>::ZeroValue();
  this->CalculateValueAndDerivative( value, derivative, false );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
{
  this->CalculateValueAndDerivative( value, derivative, true );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::CalculateValueAndDerivative( MeasureType & calculatedValue, DerivativeType & derivative, bool calculateValue ) const
{
  this->InitializeForIteration();

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->m_VirtualTransformedPointSet->GetNumberOfPoints() != this->m_FixedTransformedPointSet->GetNumberOfPoints() )
    {
    itkExceptionMacro( "Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet." );
    }

  derivative.SetSize( this->GetNumberOfParameters() );
  if( ! this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
    {
    derivative.SetSize( PointDimension * this->m_FixedTransformedPointSet->GetNumberOfPoints() );
    }
  derivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  /*
   * Split pointset in nWorkUnits ranges and sum individually
   * This splitting is required in order to avoid having the threads
   * repeatedly write to same location causing false sharing
   */
  //Use STL container to make sure no unesecarry checks are performed
  using FixedTransformedVectorContainer = typename FixedPointsContainer::STLContainerType;
  using VirtualPointsContainer = typename VirtualPointSetType::PointsContainer;
  using VirtualVectorContainer =  typename VirtualPointsContainer::STLContainerType;
  const VirtualVectorContainer &virtualTransformedPointSet =
    this->m_VirtualTransformedPointSet->GetPoints()->CastToSTLConstContainer();
  const FixedTransformedVectorContainer &fixedTransformedPointSet =
    this->m_FixedTransformedPointSet->GetPoints()->CastToSTLConstContainer();

  PointIdentifierRanges ranges = this->CreateRanges();
  std::vector< CompensatedSummation< MeasureType > > threadValues( ranges.size() );
  using CompensatedDerivative = typename std::vector< CompensatedSummation<ParametersValueType> >;
  std::vector< CompensatedDerivative > threadDerivatives( ranges.size() );
  std::function< void(unsigned int) > sumNeighborhoodValues =
      [ this, &derivative, &threadDerivatives, &threadValues, &ranges, &calculateValue,
        &virtualTransformedPointSet, &fixedTransformedPointSet]
      (unsigned int rangeIndex)
    {
    MovingTransformJacobianType  jacobian( MovingPointDimension, this->GetNumberOfLocalParameters() );
    MovingTransformJacobianType  jacobianCache;

    DerivativeType threadLocalTransformDerivative( this->GetNumberOfLocalParameters() );
    threadLocalTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

    CompensatedDerivative threadDerivativeSum( this->GetNumberOfLocalParameters() );

    CompensatedSummation< MeasureType > threadValue;
    PixelType pixel;
    NumericTraits<PixelType>::SetLength( pixel, 1 );
    for( PointIdentifier index=ranges[rangeIndex].first; index < ranges[rangeIndex].second; index++)
      {
      MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();
      LocalDerivativeType pointDerivative;

      /* Verify the virtual point is in the virtual domain.
       * If user hasn't defined a virtual space, and the active transform is not
       * a displacement field transform type, then this will always return true. */
      if( ! this->IsInsideVirtualDomain( virtualTransformedPointSet[index] ) )
        {
        continue;
        }

      if( this->m_UsePointSetData )
        {
        bool doesPointDataExist = this->m_FixedPointSet->GetPointData( index, &pixel );
        if( ! doesPointDataExist )
          {
          itkExceptionMacro( "The corresponding data for point with id " << index << " does not exist." );
          }
        }

      if( calculateValue )
        {
        this->GetLocalNeighborhoodValueAndDerivative( index, pointValue, pointDerivative, pixel );
        threadValue += pointValue;
        }
      else
        {
        pointDerivative = this->GetLocalNeighborhoodDerivative( index, pixel );
        }

      // Map into parameter space
      threadLocalTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

      if( this->m_CalculateValueAndDerivativeInTangentSpace )
        {
        for( DimensionType d = 0; d < PointDimension; ++d )
          {
          threadLocalTransformDerivative[d] += pointDerivative[d];
          }
        }
      else
        {
        this->GetMovingTransform()->
          ComputeJacobianWithRespectToParametersCachedTemporaries( virtualTransformedPointSet[index],
                                                                   jacobian,
                                                                   jacobianCache );

        for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
          {
          for( DimensionType d = 0; d < PointDimension; ++d )
            {
            threadLocalTransformDerivative[par] += jacobian(d, par) * pointDerivative[d];
            }
          }
        }
      // For local-support transforms, store the per-point result
      if( this->HasLocalSupport() || this->m_CalculateValueAndDerivativeInTangentSpace )
        {
        if( this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
          {
          this->StorePointDerivative( virtualTransformedPointSet[index], threadLocalTransformDerivative, derivative );
          }
        else
          {
          for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
            {
            derivative[this->GetNumberOfLocalParameters() * index + par] = threadLocalTransformDerivative[par];
            }
          }
        }
      for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
        {
        threadDerivativeSum[par] += threadLocalTransformDerivative[par];
        }
      }
      threadValues[rangeIndex] = threadValue;
      threadDerivatives[rangeIndex] = threadDerivativeSum;
    };

  //Sum per thread
  MultiThreaderBase::New()->ParallelizeArray( (PointIdentifier) 0,
                          (PointIdentifier) ranges.size(),
                           sumNeighborhoodValues, nullptr );

  //Sum thread results
  CompensatedSummation<MeasureType> value = 0;
  for(unsigned int i=0; i < threadValues.size(); i++)
    {
    value += threadValues[i];
    }
  MeasureType valueSum = value.GetSum();

  if( this->VerifyNumberOfValidPoints( valueSum, derivative ) )
    {
    // For global-support transforms, average the accumulated derivative result
    if( ! this->HasLocalSupport() && ! this->m_CalculateValueAndDerivativeInTangentSpace )
      {
      CompensatedDerivative localTransformDerivative( this->GetNumberOfLocalParameters() );
      for(unsigned int i=0; i< threadDerivatives.size(); i++)
        {
        for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
          {
          localTransformDerivative[par] += threadDerivatives[i][par];
          }
        }
      for( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
        {
         derivative[par] = localTransformDerivative[par].GetSum()
                                  / static_cast<DerivativeValueType>( this->m_NumberOfValidPoints );
         }
      }
    valueSum /= static_cast<MeasureType>( this->m_NumberOfValidPoints );
    }
  calculatedValue = valueSum;
  this->m_Value = valueSum;
}


template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::LocalDerivativeType
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodDerivative( const PointIdentifier point, const PixelType & pixel ) const
{
  MeasureType measure;
  LocalDerivativeType localDerivative;
  this->GetLocalNeighborhoodValueAndDerivative( point, measure, localDerivative, pixel );
  return localDerivative;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
const
typename PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PointIdentifierRanges
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::CreateRanges() const
{
  PointIdentifier nPoints = this->m_FixedTransformedPointSet->GetNumberOfPoints();
  PointIdentifier nWorkUnits = MultiThreaderBase::New()->GetNumberOfWorkUnits();
  if(nWorkUnits > nPoints || MultiThreaderBase::New()->GetMaximumNumberOfThreads() <= 1)
    {
    nWorkUnits = 1;
    }
  PointIdentifier startRange = 0;
  PointIdentifierRanges ranges;
  for(PointIdentifier p=1; p < nWorkUnits; ++p)
    {
    PointIdentifier endRange = (p * nPoints) / (double) nWorkUnits;
    ranges.push_back( PointIdentifierPair(startRange, endRange) );
    startRange = endRange;
    }
  ranges.push_back( PointIdentifierPair(startRange, nPoints) );

  return ranges;
}

/** PrintSelf */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
PointSetToPointSetMetric2v4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Fixed PointSet: " << this->m_FixedPointSet.GetPointer() << std::endl;
  os << indent << "Fixed Transform: " << this->m_FixedTransform.GetPointer() << std::endl;
  os << indent << "Moving PointSet: " << this->m_MovingPointSet.GetPointer() << std::endl;
  os << indent << "Moving Transform: " << this->m_MovingTransform.GetPointer() << std::endl;

  os << indent << "Store derivative as sparse field = ";
  if( this->GetStoreDerivativeAsSparseFieldForLocalSupportTransforms() )
    {
    os << "true." << std::endl;
    }
  else
    {
    os << "false." << std::endl;
    }

  os << indent << "Calculate in tangent space = ";
  if( this->GetCalculateValueAndDerivativeInTangentSpace() )
    {
    os << "true." << std::endl;
    }
  else
    {
    os << "false." << std::endl;
    }
}
} // end namespace itk

#endif
