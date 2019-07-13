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

#include "itkPointSetOptimalTransportMethod.h"

#include "TransportLPSolver.h"
#include "PropagationStrategy.h"
#include "NeighborhoodStrategy.h"


#include "IKMTree.h"

namespace itk
{

/** \class PointSetOptimalTransportMethod
 * \brief Base class for computing optimal transport on PointSets.
 *
 *
 * \ingroup ITKOptimalTransport
 */
template< typename TSourcePointSet, typename TTargetPointSet, typename TValue >
class ITK_TEMPLATE_EXPORT PointSetMultiscaleOptimalTransportMethod :
  public PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue > {
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(PointSetMultiscaleOptimalTransportMethod);

  /** Standard class type aliases. */
  using Self = PointSetMultiscaleOptimalTransportMethod;
  using Superclass = PointSetOptimalTransportMethod< TSourcePointSet, TTargetPointSet, TValue >;
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

  using PropagationStrategyType = PropagationStrategy<TValue>;
  using NeighborhoodStrategyType = NeighborhoodStrategy<TValue>;
  using TransportType = TransportLPSolver<double>::TransportType;

  using TransportCouplingType = typename Superclass::TransportCouplingType;

  itkSetMacro(MatchScale, bool);
  itkBooleanMacro(MatchScale);

  itkSetMacro(ScaleMass, bool);
  itkBooleanMacro(ScaleMass);

  itkSetMacro(TransportType, TransportType);
  itkGetMacro(TransportType, TransportType);

  itkSetMacro(Lambda, double);
  itkGetMacro(Lambda, double);

  itkSetMacro(MassCost, double);
  itkGetMacro(MassCost, double);

  itkSetMacro(SourceEpsilon, double);
  itkGetMacro(SourceEpsilon, double);

  itkSetMacro(TargetEpsilon, double);
  itkGetMacro(TargetEpsilon, double);

  itkSetMacro(NumberOfScalesTarget, int);
  itkGetMacro(NumberOfScalesTarget, int);

  itkSetMacro(NumberOfScalesSource, int);
  itkGetMacro(NumberOfScalesSource, int);
  
  void AddNeighborhoodPropagationStrategy(NeighborhoodStrategyType *strategy)
    {
    m_NeighborhoodStrategies.push_back(strategy);
    }

  typename TransportCouplingType::Pointer GetCoupling()
    {
    return static_cast< TransportCouplingType * >( this->ProcessObject::GetOutput(0) );
    }

protected:
  PointSetMultiscaleOptimalTransportMethod();
  virtual ~PointSetMultiscaleOptimalTransportMethod();
  void PrintSelf(std::ostream & os, Indent indent) const override;

  void GenerateData() override;

private:

  using StoppingCriterium = IKMTree<double>::StoppingCriterium;
  using SplitCriterium    = IKMTree<double>::SplitCriterium;

  /**
   * Multiscale transport solver settings
   */
  LPSolver *m_Solver;

  std::vector< NeighborhoodStrategyType* > m_NeighborhoodStrategies;
  
  PropagationStrategyType *m_PropagationStrategy1;
  PropagationStrategyType *m_PropagationStrategy2;
  
  bool m_MatchScale;
  bool m_ScaleMass;
  
  TransportType m_TransportType;
  
  double m_Lambda;
  double m_MassCost;
  double m_Exponent;
  
  int m_NumberOfScalesSource;
  int m_NumberOfScalesTarget;
  int m_MaxNeighborhoodSize;

  /**
   * GMRA (mutliscale point set representation settings)
   */
  SplitCriterium    m_SourceSplitCriterium;
  StoppingCriterium m_SourceStoppingCriterium;
  
  int m_SourceMinimumPoints;
  int m_SourceMaxIterations;
  int m_SourceNumberOfKids;
  
  double m_SourceEpsilon;
  double m_SourceThreshold;

  SplitCriterium    m_TargetSplitCriterium;
  StoppingCriterium m_TargetStoppingCriterium;
  
  int m_TargetMinimumPoints;
  int m_TargetMaxIterations;
  int m_TargetNumberOfKids;
  
  double m_TargetEpsilon;
  double m_TargetThreshold;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetMultiscaleOptimalTransportMethod.hxx"
#endif

#endif
