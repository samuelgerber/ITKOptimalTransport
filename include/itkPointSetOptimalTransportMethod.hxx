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

template< typename TSourcePointSet, typename TTargetPointSet >
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::PointSetOptimalTransportMethod()
{
  this->SetNumberOfRequiredOutputs(1);


  TransportPlanOutputPointer transformDecorator =
    itkDynamicCastInDebugMode< TransportPlanType * >(this->MakeOutput(0).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );
}


template< typename TSourcePointSet, typename TTargetPointSet >
void
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
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
  auto * transformOutput = static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );

  transformOutput->Set( m_Transform );
}


template< typename TSourcePointSet, typename TTargetPointSet >
const typename PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >::TransformOutputType *
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::GetOutput() const
{
  return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
}

template< typename TSourcePointSet, typename TTargetPointSet >
DataObject::Pointer
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::MakeOutput(DataObjectPointerArraySizeType output)
{
  if (output > 0)
  {
    itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs.");
  }
  return TransformOutputType::New().GetPointer();
}

template< typename TSourcePointSet, typename TTargetPointSet >
ModifiedTimeType
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
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

template< typename TSourcePointSet, typename TTargetPointSet >
void
PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  itkPrintSelfObjectMacro( SourcePointSet );
  itkPrintSelfObjectMacro( TargetPointSet );

}
} // end namespace itk
#endif
