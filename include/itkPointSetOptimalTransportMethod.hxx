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
#ifndef itkPointSetOptimalTransportMethod_hxx
#define itkPointSetOptimalTransportMethod_hxx

#include "itkPointSetOptimalTransportMethod.h"

namespace itk
{

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PointSetOptimalTransportMethod()
{
  this->SetNumberOfRequiredOutputs(1);


  TransportPlanOutputPointer output =
    itkDynamicCastInDebugMode< TransportPlanType * >(this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}


template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::Initialize()
{
  if ( !m_SourcePointSet )
    {
    itkExceptionMacro(<< "SourcePointSet is not present");
    }

  if ( !m_TargetPointSet )
    {
    itkExceptionMacro(<< "TargetPointSet is not present");
    }


  // Connect the transform to the Decorator
  auto * transportOutput = static_cast< TransportPlanType * >( this->ProcessObject::GetOutput(0) );

  transportOutput->Set( m_TransportPlan );
}


template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
const typename PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >::TransformOutputType *
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::GetOutput() const
{
  return static_cast< const TransportPlanType * >( this->ProcessObject::GetOutput(0) );
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
DataObject::Pointer
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::MakeOutput(DataObjectPointerArraySizeType output)
{
  if (output > 0)
  {
    itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs.");
  }
  return TransportPlanType::New().GetPointer();
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
ModifiedTimeType
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::GetMTime() const
{
  ModifiedTimeType mtime = Superclass::GetMTime();
  ModifiedTimeType m;

  // Some of the following should be removed once ivars are put in the
  // input and output lists
  if ( m_SourcePointSet )
    {
    m = m_SourcePointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_TargetPointSet )
    {
    m = m_TargetPointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  return mtime;
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  itkPrintSelfObjectMacro( SourcePointSet );
  itkPrintSelfObjectMacro( TargetPointSet );

}
} // end namespace itk
#endif
