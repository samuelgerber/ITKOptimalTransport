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
#include "PointSetGMRADataObject.h"
#include "IKMTree.h"

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
  //Create Source GMRA object
  PointSetGMRADataObject source(m_SourcePointSet);
  std::vector<int> sourcePts(source.numberOfPoints() );
  for(unsigned int i=0; i<sourcePts.size(); i++){
    sourcePts[i] = i;
  };
  IKMTree<double> *gmraSource = new IKMTree<double>(source);
  gmra->setStoppingCriterium( m_SourceStoppingCriterium );
  gmra->setSplitCriterium( m_SourceSplitCriterium );
  gmra->dataFactory = new L2GMRAKmeansDataFactory<double>();
  gmra->epsilon = m_SourceEpsilon;
  gmra->nKids = m_SourceNKids;
  gmra->threshold = m_SourceThreshold;
  gmra->maxIter = m_SourceMaxIterations;
  gmra->minPoints = m_SourceMinimumPoints;
  gmra->addPoints(sourcePts);


  //Create Target GMRA object
  PointSetGMRADataObject target(m_TargetPointSet);
  std::vector<int> targetPts(target.numberOfPoints() );
  for(unsigned int i=0; i<targetPts.size(); i++){
    targetPts[i] = i;
  };
  IKMTree<double> *gmraTarget = new IKMTree<double>(target);
  gmra->setStoppingCriterium( m_TargetStoppingCriterium );
  gmra->setSplitCriterium( m_TargetSplitCriterium );
  gmra->dataFactory = new L2GMRAKmeansDataFactory<double>();
  gmra->epsilon = m_TargetEpsilon;
  gmra->nKids = m_TargetNKids;
  gmra->threshold = m_TargetThreshold;
  gmra->maxIter = m_TargetMaxIterations;
  gmra->minPoints = m_TargetMinimumPoints;
  gmra->addPoints(targetPts);

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
