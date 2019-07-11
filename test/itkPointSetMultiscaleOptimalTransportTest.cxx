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
#include "itkPointSetMultiscaleOptimalTransportMethod.h"
#include "itkOptimalTransportPointSetMetric.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkCommand.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"
#include <fstream>
#include <iostream>

template<typename TFilter>
class itkPointSetMetricRegistrationTestCommandIterationUpdate : public itk::Command
{
public:
  using Self = itkPointSetMetricRegistrationTestCommandIterationUpdate;

  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro( Self );

protected:
  itkPointSetMetricRegistrationTestCommandIterationUpdate() = default;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *) caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }
    const auto * optimizer = dynamic_cast< const TFilter * >( object );

    if( !optimizer )
      {
      itkGenericExceptionMacro( "Error dynamic_cast failed" );
      }
    std::cout << "It: " << optimizer->GetCurrentIteration() << " metric value: " << optimizer->GetCurrentMetricValue();
    std::cout << std::endl;
    }
};



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


  using PixelType = double;
  using FixedImageType = itk::Image<PixelType, Dimension>;
  using MovingImageType = itk::Image<PixelType, Dimension>;

  FixedImageType::SizeType fixedImageSize;
  FixedImageType::PointType fixedImageOrigin;
  FixedImageType::DirectionType fixedImageDirection;
  FixedImageType::SpacingType fixedImageSpacing;

  fixedImageSize.Fill( 600 );
  fixedImageOrigin.Fill( -300 );
  fixedImageDirection.SetIdentity();
  fixedImageSpacing.Fill( 1 );

  FixedImageType::Pointer fixedImage = FixedImageType::New();
  fixedImage->SetRegions( fixedImageSize );
  fixedImage->SetOrigin( fixedImageOrigin );
  fixedImage->SetDirection( fixedImageDirection );
  fixedImage->SetSpacing( fixedImageSpacing );
  fixedImage->Allocate();

  using AffineTransformType = itk::AffineTransform<double, Dimension>;
  AffineTransformType::Pointer transform = AffineTransformType::New();
  transform->SetIdentity();

  using PointSetMetricType = itk::OptimalTransportPointSetMetric<PointSetType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetCoupling( coupling );
  metric->SetMovingTransform( transform );
  metric->SetFixedPointSet( fixedPoints );
  metric->SetMovingPointSet( movingPoints );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->Initialize();

  // scales estimator
  using RegistrationParameterScalesFromShiftType =
          itk::RegistrationParameterScalesFromPhysicalShift< PointSetMetricType >;
  RegistrationParameterScalesFromShiftType::Pointer shiftScaleEstimator =
          RegistrationParameterScalesFromShiftType::New();
  shiftScaleEstimator->SetMetric( metric );
  // needed with pointset metrics
  shiftScaleEstimator->SetVirtualDomainPointSet( metric->GetVirtualTransformedPointSet() );

  // optimizer
  using OptimizerType = itk::GradientDescentOptimizerv4;
  OptimizerType::Pointer  optimizer = OptimizerType::New();
  optimizer->SetMetric( metric );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetScalesEstimator( shiftScaleEstimator );
  optimizer->SetMaximumStepSizeInPhysicalUnits( 3 );

  using CommandType = itkPointSetMetricRegistrationTestCommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  optimizer->SetMinimumConvergenceValue( 0.0 );
  optimizer->SetConvergenceWindowSize( 10 );

  using AffineRegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, AffineTransformType, FixedImageType, PointSetType>;
  AffineRegistrationType::Pointer affineSimple = AffineRegistrationType::New();
  affineSimple->SetNumberOfLevels(1);
  affineSimple->SetObjectName( "affineSimple" );
  affineSimple->SetFixedPointSet( fixedPoints );
  affineSimple->SetMovingPointSet( movingPoints );
  affineSimple->SetInitialTransform( transform );
  affineSimple->SetMetric( metric );
  affineSimple->SetOptimizer( optimizer );

   try
    {
    std::cout << "Trimmed point set affine registration update" << std::endl;
    affineSimple->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }
/*
  if(metric->GetValue() < 2){
    return EXIT_SUCCESS;
  }
  else{
    return EXIT_FAILURE;
  }
*/

  return EXIT_SUCCESS;
}
