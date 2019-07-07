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
#ifndef itkPointSetOptimalTransportMethod_h
#define itkPointSetOptimalTransportMethod_h

#include "itkProcessObject.h"
#include "itkPointSet.h"
#include "itkDataObject.h"
#include "itkTransportCoupling.h"


namespace itk
{

/** \class PointSetOptimalTransportMethod
 * \brief Base class for computing optimal transport on PointSets.
 *
 *
 * \ingroup ITKOptimalTransport
 */
template< typename TSourcePointSet, typename TTargetPointSet, typename TValue = double >
class ITK_TEMPLATE_EXPORT PointSetOptimalTransportMethod : public ProcessObject
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(PointSetOptimalTransportMethod);

  /** Transport Map Definition **/
  using SourceCoordRepType = typename TSourcePointSet::CoordRepType;
  using SourcePointIdentifier = typename TSourcePointSet::PointIdentifier;
  using SourcetPointType = typename TSourcePointSet::PointType;

  using TargetCoordRepType = typename TTargetPointSet::CoordRepType;
  using TargetPointIdentifier = typename TTargetPointSet::PointIdentifier;
  using TargetPointType = typename TTargetPointSet::PointType;
  
  using TransportCouplingType = TransportCoupling<SourcePointIdentifier, TargetPointIdentifier, TValue>;
  using TransportCouplingPointer = typename TransportCouplingType::Pointer;

  /** Standard class type aliases. */
  using Self = PointSetOptimalTransportMethod;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetOptimalTransportMethod, ProcessObject);

  /**  Type of the Fixed PointSet. */
  using SourcePointSetType = TSourcePointSet;
  using SourcePointSetConstPointer = typename SourcePointSetType::ConstPointer;

  /**  Type of the Moving PointSet. */
  using TargetPointSetType = TTargetPointSet;
  using TargetPointSetConstPointer = typename TargetPointSetType::ConstPointer;

  /** Smart Pointer type to a DataObject. */
  using DataObjectPointer = typename DataObject::Pointer;

  /** Set/Get the Fixed PointSet. */
  itkSetConstObjectMacro(SourcePointSet, SourcePointSetType);
  itkGetConstObjectMacro(SourcePointSet, SourcePointSetType);

  /** Set/Get the Moving PointSet. */
  itkSetConstObjectMacro(TargetPointSet, TargetPointSetType);
  itkGetConstObjectMacro(TargetPointSet, TargetPointSetType);




  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) override;

  ModifiedTimeType GetMTime() const override;

protected:
  PointSetOptimalTransportMethod();
  ~PointSetOptimalTransportMethod() override = default;
  void PrintSelf(std::ostream & os, Indent indent) const override;

  virtual void Initialize();
  void GenerateData() override;

private:
  TargetPointSetConstPointer m_TargetPointSet;
  SourcePointSetConstPointer  m_SourcePointSet;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetOptimalTransportMethod.hxx"
#endif

#endif

