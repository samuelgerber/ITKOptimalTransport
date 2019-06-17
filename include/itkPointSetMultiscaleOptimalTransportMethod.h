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
#ifndef itkPointSetMultiscaleOptimalTransportMethod_h
#define itkPointSetMultiscaleOptimalTransportMethod_h

#include "itkPointSetMultiscaleOptimalTransportMethod.h"

namespace itk

/** \class PointSetOptimalTransportMethod
 * \brief Base class for computing optimal transport on PointSets.
 *
 *
 * \ingroup ITKOptimalTransport
 */
template< typename TSourcePointSet, typename TTargetPointSet >
class ITK_TEMPLATE_EXPORT PointSetMultiscaleOptimalTransportMethod : 
  public PointSetOptimalTransportMethod<TSourcePointSet, TTargetPointSet> {
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(PointSetMultiscaleOptimalTransportMethod);

  /** Standard class type aliases. */
  using Self = PointSetMultiscaleOptimalTransportMethod;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetMultiscaleOptimalTransportMethod, ProcessObject);

  /**  Type of the Fixed PointSet. */
  using SourcePointSetType = TSourcePointSet;
  using SourcePointSetConstPointer = typename SourcePointSetType::ConstPointer;

  /**  Type of the Moving PointSet. */
  using TargetPointSetType = TTargetPointSet;
  using TargetPointSetConstPointer = typename TargetPointSetType::ConstPointer;


protected:
  PointSetMultiscaleOptimalTransportMethod();
  ~PointSetMultiscaleOptimalTransportMethod() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  void GenerateData() override;


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetMultiscaleOptimalTransportMethod.hxx"
#endif

#endif

