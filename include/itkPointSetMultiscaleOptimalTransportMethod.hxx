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

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PointSetMultiscaleOptimalTransportMethod()
{
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::Initialize()
{
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::GenerateData()
{
}



template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
ModifiedTimeType
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::GetMTime() const
{
  ModifiedTimeType mtime = Superclass::GetMTime();
  ModifiedTimeType m;

  return mtime;
}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk
#endif
