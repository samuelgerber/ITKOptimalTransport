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

#include "itkMersenneTwisterRandomVariateGenerator.h"
//#include "itkPointSetOptimalTransportMethod.h"
#include "../include/itkPointSetMultiscaleOptimalTransportMethod.h"
//#include "itkPointSetMultiscaleOptimalTransportMethod.h"

#include <fstream>
#include <iostream>




int itkPointSetMultiscaleOptimalTransportTest( int argc, char *argv[] )
{
  constexpr unsigned int Dimension = 2;

  unsigned int numberOfIterations = 100;
  if( argc > 1 )
    {
    numberOfIterations = std::stoi( argv[1] );
    }

  using PointSetType = itk::PointSet<double, Dimension>;

  using PointType = PointSetType::PointType;

  PointSetType::Pointer fixedPoints = PointSetType::New();
  fixedPoints->Initialize();

  PointSetType::Pointer movingPoints = PointSetType::New();
  movingPoints->Initialize();

  //itk::MultiThreaderBase::New()->SetMaximumNumberOfThreads( 8 );
  //itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads( 8 );
  std::cout << "MaxNumberOfThreads: ";
  std::cout << itk::MultiThreaderBase::New()->GetMaximumNumberOfThreads() << std::endl;
  std::cout << "NumberOfWorkUnits: ";
  std::cout << itk::MultiThreaderBase::New()->GetNumberOfWorkUnits() << std::endl;

  using GeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize(2019);

  // Generate two noisy ellipses
  unsigned int nSourcePoints= 1000;
  for(int i=0; i< nSourcePoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 2*itk::Math::pi);
    PointType fixedPoint;
    fixedPoint[0] = 0.5 * radius * std::cos( theta ) + generator->GetNormalVariate();
    fixedPoint[1] = 2 * radius * std::sin( theta ) + generator->GetNormalVariate();
    fixedPoints->SetPoint( i, fixedPoint );
    }

  unsigned int nTargetPoints= 1200;
  for(int i=0; i< nTargetPoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 1*itk::Math::pi);
    PointType movingPoint;
    movingPoint[0] = 0.75 * radius * std::cos( theta ) + generator->GetNormalVariate();
    movingPoint[1] = 1.5 * radius * std::sin( theta ) + generator->GetNormalVariate();
    movingPoints->SetPoint( i, movingPoint );
    }

  using OptimalTransportType = itk::PointSetMultiscaleOptimalTransportMethod<PointSetType, PointSetType, double>;

  OptimalTransportType::Pointer ot = OptimalTransportType::New();
  ot->SetSourcePointSet( fixedPoints );
  ot->SetTargetPointSet( movingPoints );

  ot->Update();

  std::cout << "Update called" << std::endl;

  typename OptimalTransportType::TransportCouplingType::Pointer coupling = ot->GetCoupling();
  std::cout << coupling << std::endl;

  return EXIT_SUCCESS;
}
