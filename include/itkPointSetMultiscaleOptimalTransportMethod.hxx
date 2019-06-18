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
#include "LemonSolver.h"

namespace itk
{

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PointSetMultiscaleOptimalTransportMethod()
{
  m_Solver = new LemonSolver();

}

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::~PointSetMultiscaleOptimalTransportMethod()
{
  delete m_Solver;

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
  std::vector<double> sourceWeights(source.numberOfPoints() );
  for(unsigned int i=0; i<sourcePts.size(); i++){
    sourcePts[i] = i;
    sourceWeights[i] = 1;
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
  std::vector<double> targetWeights(target.numberOfPoints() );
  for(unsigned int i=0; i<targetPts.size(); i++){
    targetPts[i] = i;
    targetWeights[i] = 1.0;
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

  
  NodeDistance<double> *dist = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
  gmraSource->computeRadii(dist);
  gmraSource->computeLocalRadii(dist);
  
  gmraTarget->computeRadii(dist);
  gmraTarget->computeLocalRadii(dist);

  GenericGMRANeighborhood<double> sourceNeighborhood(gmraSource, dist);
  GenericGMRANeighborhood<double> targetNeighborhood(gmraTarget, dist);

  std::vector< MultiscaleTransportLevel<double> * > sourceLevels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(sourceNeighborhood, sourceWeights, false);

  std::vector< MultiscaleTransportLevel<double> * > targetLevels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(targetNeighborhood, targetWeights, false);

  TransportLPSolver<double> *trpSolver =
      new TransportLPSolver<double>( m_Solver, m_TransportType, m_MassCost, m_Lambda );
  MultiscaleTransportLP<double> transport = new MultiscaleTransportLP<double>(&trpSolver);
  transport.setPropagationStrategy1(m_PropagationStrategy1);
  transport.setPropagationStrategy1(m_PropagationStrategy2);
  for(int i=0; i< m_NeighborhoodPropagations.size(); i++){
    transport.addNeighborhoodPropagationStrategy( m_NeighborhoodPropagations[i] );
  }

  std::vector< TransportPlan<double> * > sols = transport.solve( sourceLevels, targetLevels, 
      m_Exponent, m_NumberOfScalesSource, m_NumberOfScalesTarget, m_MatchScale, m_ScaleMass);

  delete gmraSource;
  delete gmraTarget;
  delete dist; 

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
