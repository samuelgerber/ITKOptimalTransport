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
#ifndef itkPointSetMultiscaleOptimalTransportMethod_hxx
#define itkPointSetMultiscaleOptimalTransportMethod_hxx

#include "itkPointSetMultiscaleOptimalTransportMethod.h"

namespace itk
{

template< typename TSourcePointSet, typename TTargetPointSet >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::PointSetMultiscaleOptimalTransportMethod()
{
;
}

template< typename TSourcePointSet, typename TTargetPointSet >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::Initialize()
{

}

template< typename TSourcePointSet, typename TTargetPointSet >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::GenerateData()
{

}

template< typename TSourcePointSet, typename TTargetPointSet >
const typename PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >::TransformOutputType *
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::GetOutput() const
{
  return static_cast< const TransportPlan * >( this->ProcessObject::GetOutput(0) );
}

template< typename TSourcePointSet, typename TTargetPointSet >
DataObject::Pointer
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::MakeOutput(DataObjectPointerArraySizeType output)
{
  if (output > 0)
  {
    itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs.");
  }
  return TransportPlan::New().GetPointer();
}

template< typename TSourcePointSet, typename TTargetPointSet >
ModifiedTimeType
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::GetMTime() const
{
  ModifiedTimeType mtime = Superclass::GetMTime();
  ModifiedTimeType m;

  return mtime;
}

template< typename TSourcePointSet, typename TTargetPointSet >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

}
} // end namespace itk
#endif
