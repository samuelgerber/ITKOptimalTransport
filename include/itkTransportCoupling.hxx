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
#ifndef itkTransportCoupling_hxx
#define itkTransportCoupling_hxx

#include "itkTransportCoupling.h"
/*
 *\ingroup OptimalTransport"
 */
namespace itk
{
template<typename TSourcePointIdentifier, typename TTargetPointIdentifier, typename TValue>
TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>
::TransportCoupling()
{

}

template<typename TSourcePointIdentifier, typename TTargetPointIdentifier, typename TValue>
TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>
::~TransportCoupling()
{

}

template<typename TSourcePointIdentifier, typename TTargetPointIdentifier, typename TValue>
void
TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "TransportCoupling: ";
  os << indent << m_Map.size() << " sources";
}


template<typename TSourcePointIdentifier, typename TTargetPointIdentifier, typename TValue>
void
TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>
::AddPath(TSourcePointIdentifier source, TTargetPointIdentifier target, TValue weight)
{
  m_Map[source][target] = weight;
}


template<typename TSourcePointIdentifier, typename TTargetPointIdentifier, typename TValue>
typename TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>::TransportMap &
TransportCoupling<TSourcePointIdentifier, TTargetPointIdentifier, TValue>
::GetMap()
{
  return m_Map;
}

} // end of namespace itk

#endif
