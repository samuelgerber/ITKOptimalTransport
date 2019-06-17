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
#ifndef itkTransportPlan_h
#define itkTransportPlan_h

#include <queue>
#include <vector>

#include "itkPoint.h"
#include "itkSize.h"
#include "itkObject.h"
#include "itkArray.h"



namespace itk
{

/** \class TransportPlam
 *  \brief
 *
 *  */

template<typename TSourePointIdentifier, typename TTargetPointIdentifier, typename TValue = double>
class ITK_TEMPLATE_EXPORT TransportPlan:public Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(TransportPlan);

  /** Standard class type aliases */
  using Self = TransportPlan;
  using Superclass = Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information (and related methods) */
  itkTypeMacro(TransportPlan, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  using std::pair< TSourePointIdentifier, TTargetPointIdentifier> TransportPath;
  using std::map< TransportPath, TValue> TransportMap;

  void AddPath(TSourePointIdentifier source, TTargetPointIdentifier target, TValue weight);

  TransportMap &GetMap();

protected:
  /** Constructor */
  TransportPlan();
  ~TransportPlan() override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;



private:

  TransportMap m_Map;

};  // end of class

} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTransportPlan.hxx"
#endif

#endif
