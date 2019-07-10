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
#include "IteratedCapacityPropagationStrategy.h"
#include "EigenEuclideanMetric.h"
#include "GMRANeighborhood.h"
#include "GMRAMultiscaleTransport.h"
#include "MultiscaleTransportLP.h"

namespace itk
{

template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::PointSetMultiscaleOptimalTransportMethod()
{
  m_Solver = new LemonSolver();
  m_PropagationStrategy1 = new IteratedCapacityPropagationStrategy<TValue>(3, 0);
  m_PropagationStrategy2 = new IteratedCapacityPropagationStrategy<TValue>(3, 0);
  m_MaxNeighborhoodSize = NumericTraits<int>::max();
  m_MatchScale = false;
  m_ScaleMass = true;
  m_Lambda = 0;
  m_MassCost = 0;
  m_Exponent = 2;
  m_NumberOfScalesSource = -1;
  m_NumberOfScalesTarget= -1;
  m_TransportType = TransportLPSolver<double>::BALANCED;

  m_SourceSplitCriterium = IKMTree<TValue>::ADAPTIVE_FIXED;
  m_SourceStoppingCriterium = IKMTree<TValue>::RELATIVE_RADIUS;
  //TODO need to fix epsilon set to 0 bug
  m_SourceEpsilon = 0.00000001;
  m_SourceNumberOfKids = 8;
  m_SourceThreshold = 0;
  m_SourceMaxIterations = 100;
  m_SourceMinimumPoints = 1;

  m_TargetSplitCriterium = IKMTree<TValue>::ADAPTIVE_FIXED;
  m_TargetStoppingCriterium = IKMTree<TValue>::RELATIVE_RADIUS;
  m_TargetEpsilon = 0.00000001;
  m_TargetNumberOfKids = 8;
  m_TargetThreshold = 0;
  m_TargetMaxIterations = 100;
  m_TargetMinimumPoints = 1;
}


template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::~PointSetMultiscaleOptimalTransportMethod()
{
  delete m_Solver;
  delete m_PropagationStrategy1;
  delete m_PropagationStrategy2;
}


template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
void
PointSetMultiscaleOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >
::GenerateData()
{

  //Create Source GMRA object
  PointSetGMRADataObject<TSourcePointSet> source( this->GetSourcePointSet() );
  std::vector<int> sourcePts( source.numberOfPoints() );
  std::vector<double> sourceWeights( source.numberOfPoints() );
  for(unsigned int i=0; i<sourcePts.size(); i++){
    sourcePts[i] = i;
    sourceWeights[i] = 1.0;
  };

  IKMTree<double> *gmraSource = new IKMTree<double>( &source );
  gmraSource->setStoppingCriterium( m_SourceStoppingCriterium );
  gmraSource->setSplitCriterium( m_SourceSplitCriterium );
  gmraSource->dataFactory = new L2GMRAKmeansDataFactory<double>();
  gmraSource->epsilon = m_SourceEpsilon;
  gmraSource->nKids = m_SourceNumberOfKids;
  gmraSource->threshold = m_SourceThreshold;
  gmraSource->maxIter = m_SourceMaxIterations;
  gmraSource->minPoints = m_SourceMinimumPoints;
  gmraSource->addPoints( sourcePts );


  //Create Target GMRA object
  PointSetGMRADataObject<TTargetPointSet> target( this->GetTargetPointSet() );
  std::vector<int> targetPts( target.numberOfPoints() );
  std::vector<double> targetWeights( target.numberOfPoints() );
  for(unsigned int i=0; i<targetPts.size(); i++){
    targetPts[i] = i;
    targetWeights[i] = 1.0;
  };
  IKMTree<double> *gmraTarget = new IKMTree<double>( &target );
  gmraTarget->setStoppingCriterium( m_TargetStoppingCriterium );
  gmraTarget->setSplitCriterium( m_TargetSplitCriterium );
  gmraTarget->dataFactory = new L2GMRAKmeansDataFactory<double>();
  gmraTarget->epsilon = m_TargetEpsilon;
  gmraTarget->nKids = m_TargetNumberOfKids;
  gmraTarget->threshold = m_TargetThreshold;
  gmraTarget->maxIter = m_TargetMaxIterations;
  gmraTarget->minPoints = m_TargetMinimumPoints;
  gmraTarget->addPoints( targetPts );

  std::cout << targetPts.size() << std::endl;

  NodeDistance<double> *distS = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
  gmraSource->computeRadii(distS);
  gmraSource->computeLocalRadii(distS);

  NodeDistance<double> *distT = new CenterNodeDistance<double>( new EuclideanMetric<double>() );
  gmraTarget->computeRadii(distT);
  gmraTarget->computeLocalRadii(distT);

  GenericGMRANeighborhood<double> sourceNeighborhood(gmraSource, distS);
  GenericGMRANeighborhood<double> targetNeighborhood(gmraTarget, distT);

  std::vector< MultiscaleTransportLevel<double> * > sourceLevels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(sourceNeighborhood, sourceWeights, false);

  std::vector< MultiscaleTransportLevel<double> * > targetLevels =
      GMRAMultiscaleTransportLevel<double>::buildTransportLevels(targetNeighborhood, targetWeights, false);

  std::cout << targetLevels.size() << std::endl;
  TransportLPSolver<double> *trpSolver =
          new TransportLPSolver<double>( m_Solver, m_TransportType, m_MassCost, m_Lambda );
  MultiscaleTransportLP<double> transport( trpSolver );
  transport.setPropagationStrategy1(m_PropagationStrategy1);
  transport.setPropagationStrategy1(m_PropagationStrategy2);
  for(int i=0; i< m_NeighborhoodStrategies.size(); i++)
    {
    transport.addNeighborhodStrategy( m_NeighborhoodStrategies[i] );
    }

  std::vector< TransportPlan<double> * > sols = transport.solve( sourceLevels, targetLevels,
      m_Exponent, m_NumberOfScalesSource, m_NumberOfScalesTarget, m_MatchScale, m_ScaleMass);



  auto * transportOutput = static_cast< TransportCouplingType * >( this->ProcessObject::GetOutput(0) );

  using Path = TransportPlan<double>::Path;
  TransportPlan<double> *sol = sols[sols.size()-1];
  std::cout << sols.size() << std::endl;
  for(sol->pathIteratorBegin(); ! sol->pathIteratorIsAtEnd(); sol->pathIteratorNext() )
    {
    Path &path = sol->pathIteratorCurrent();
    GMRATransportNode<double> *from = (GMRATransportNode<double> *) path.from;
    GMRATransportNode<double> *to = (GMRATransportNode<double> *) path.to;
    if(path.w > 0)
      {
      std::cout << path.w << ": " << from->getID() << " -> " << to->getID() << std::endl;
      transportOutput->AddPath(from->getID(), to->getID(), path.w);
      }
    }

  delete gmraSource;
  delete gmraTarget;
  delete distS;
  delete distT;

  for(int i = 0; i<sols.size(); i++)
    {
    delete sols[i];
    }

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
