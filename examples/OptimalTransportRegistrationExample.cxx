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
#include "itkPointSetMultiscaleOptimalTransportMethod.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkCommand.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkAffineTransform.h"
#include "itkOptimalTransportPointSetMetric.h"
#include "itkTimeProbe.h"

#include <fstream>
#include <iostream>


constexpr unsigned int Dimension = 2;
using PointSetType = itk::PointSet<double, Dimension>;
using PointType = PointSetType::PointType;

using PixelType = double;
using FixedImageType = itk::Image<PixelType, Dimension>;
using MovingImageType = itk::Image<PixelType, Dimension>;

FixedImageType::SizeType fixedImageSize;
FixedImageType::PointType fixedImageOrigin;
FixedImageType::DirectionType fixedImageDirection;
FixedImageType::SpacingType fixedImageSpacing;


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




template<typename TPointSetMetricType>
void runRegistration( PointSetType::Pointer fixedPoints,
                      PointSetType::Pointer movingPoints,
                      TPointSetMetricType *metric,
                      FixedImageType::Pointer fixedImage,
                      std::string filename){

  unsigned int numberOfIterations = 100;

  using AffineTransformType = itk::AffineTransform<double, Dimension>;
  AffineTransformType::Pointer transform = AffineTransformType::New();
  transform->SetIdentity();

  // scales estimator
  using RegistrationParameterScalesFromShiftType =
          typename itk::RegistrationParameterScalesFromPhysicalShift< TPointSetMetricType >;
  typename RegistrationParameterScalesFromShiftType::Pointer shiftScaleEstimator =
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
  optimizer->SetMaximumStepSizeInPhysicalUnits( 1 );
  optimizer->SetMinimumConvergenceValue( 0.0 );
  optimizer->SetConvergenceWindowSize( 2000 );

  using CommandType = itkPointSetMetricRegistrationTestCommandIterationUpdate<OptimizerType>;
  CommandType::Pointer observer = CommandType::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  using AffineRegistrationType = itk::ImageRegistrationMethodv4<FixedImageType,
        MovingImageType, AffineTransformType, FixedImageType, PointSetType>;
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
    //std::cout << "Trimmed point set affine registration update" << std::endl;
    affineSimple->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    }

  //Store point sets in CSV for visual inspection
  AffineTransformType::InverseTransformBasePointer affineInverseTransform =
          affineSimple->GetModifiableTransform()->GetInverseTransform();
  std::ofstream myfile;
  myfile.open (filename);
  for( unsigned int n=0; n < movingPoints->GetNumberOfPoints(); n++ )
    {
    // compare the points in virtual domain
    PointType movingPoint = movingPoints->GetPoint( n );
    PointType transformedMovingPoint =
            affineInverseTransform->TransformPoint( movingPoints->GetPoint( n ) );
    myfile << "Moving";
    for(int i=0; i<PointType::PointDimension; i++)
      {
      myfile << ", " << movingPoint[i];
      }
    myfile << std::endl;

    myfile << "MovingTransformed";
    for(int i=0; i<PointType::PointDimension; i++)
      {
      myfile << ", " << transformedMovingPoint[i];
      }
    myfile << std::endl;
    }
  for( unsigned int n=0; n < fixedPoints->GetNumberOfPoints(); n++ )
    {
    // compare the points in virtual domain
    PointType fixedPoint = fixedPoints->GetPoint( n );
    PointType transformedFixedPoint =
            affineSimple->GetModifiableTransform()->TransformPoint( fixedPoints->GetPoint( n ) );
    myfile << "Fixed";
    for(int i=0; i<PointType::PointDimension; i++)
      {
      myfile << ", " << fixedPoint[i];
      }
    myfile << std::endl;

    myfile << "FixedTransformed";
    for(int i=0; i<PointType::PointDimension; i++)
      {
      myfile << ", " << transformedFixedPoint[i];
      }
    myfile << std::endl;
    }
  myfile.close();

}

int main( int argc, char *argv[] )
{
  constexpr unsigned int Dimension = 2;

  unsigned int numberOfIterations = 500;
  if( argc > 1 )
    {
    numberOfIterations = std::stoi( argv[1] );
    }

  std::cout << "MaxNumberOfThreads: ";
  std::cout << itk::MultiThreaderBase::New()->GetMaximumNumberOfThreads() << std::endl;
  std::cout << "NumberOfWorkUnits: ";
  std::cout << itk::MultiThreaderBase::New()->GetNumberOfWorkUnits() << std::endl;

  PointSetType::Pointer fixedPoints = PointSetType::New();
  fixedPoints->Initialize();

  PointSetType::Pointer movingPoints = PointSetType::New();
  movingPoints->Initialize();


  using GeneratorType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
  GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();

  // Generate two noisy ellipses
  unsigned int nSourcePoints= 4000;
  for(int i=0; i< nSourcePoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 2*itk::Math::pi);
    PointType fixedPoint;
    fixedPoint[0] = 0.5 * radius * std::cos( theta ) + generator->GetNormalVariate();
    fixedPoint[1] = 2 * radius * std::sin( theta ) + generator->GetNormalVariate();
    fixedPoints->SetPoint( i, fixedPoint );
    }

  unsigned int nTargetPoints= 401;
  for(int i=0; i< nTargetPoints; i++ )
    {
    float radius = 100.0;

    float theta = generator->GetUniformVariate(0, 1*itk::Math::pi);
    PointType movingPoint;
    movingPoint[0] = 0.75 * radius * std::cos( theta ) + generator->GetNormalVariate();
    movingPoint[1] = 1.5 * radius * std::sin( theta ) + generator->GetNormalVariate();
    movingPoints->SetPoint( i, movingPoint );
    }

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


  std::cout << "Running examples with " << nSourcePoints << " source points and " << nTargetPoints << " target points." << std::endl;
  std::cout << "Using a maximum of " << itk::MultiThreaderBase::New()->GetMaximumNumberOfThreads() << " threads and ";
  std::cout << itk::MultiThreaderBase::New()->GetNumberOfWorkUnits() << " work units." << std::endl;

  {
  itk::TimeProbe clock;
  clock.Start();

  using OptimalTransportType = itk::PointSetMultiscaleOptimalTransportMethod<PointSetType, PointSetType, double>;
  OptimalTransportType::Pointer ot = OptimalTransportType::New();
  ot->SetSourcePointSet( fixedPoints );
  ot->SetTargetPointSet( movingPoints );
  ot->Update();
  typename OptimalTransportType::TransportCouplingType::Pointer coupling = ot->GetCoupling();
  coupling->SaveToCsv( "ot-balanced-coupling.csv" );

  using PointSetMetricType = itk::OptimalTransportPointSetMetric<PointSetType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetMovingTransform( transform );
  metric->SetFixedPointSet( fixedPoints );
  metric->SetMovingPointSet( movingPoints );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->SetCoupling( coupling );
  metric->Initialize();
  runRegistration<PointSetMetricType>(fixedPoints, movingPoints, metric, fixedImage, "ot-balanced.csv");
  clock.Stop();
  std::cout << "OT1 Metric Time : " << clock.GetTotal() << std::endl;
  }


  {
  itk::TimeProbe clock;
  clock.Start();

  PointSetType::Pointer fixedTransformedPoints = fixedPoints;
  for(int i=0; i<3; i++){
  using OptimalTransportType = itk::PointSetMultiscaleOptimalTransportMethod<PointSetType, PointSetType, double>;
  OptimalTransportType::Pointer ot = OptimalTransportType::New();
  ot->SetSourcePointSet( fixedTransformedPoints );
  ot->SetTargetPointSet( movingPoints );
  ot->SetScaleMass(false);
  ot->SetTransportType( TransportLPSolver<double>::UNBALANCED_FREE );
  ot->SetMassCost( 10000 );
  ot->Update();
  typename OptimalTransportType::TransportCouplingType::Pointer coupling = ot->GetCoupling();
  coupling->SaveToCsv( "ot-unbalanced-coupling.csv" );

  using PointSetMetricType = itk::OptimalTransportPointSetMetric<PointSetType>;
  PointSetMetricType::Pointer metric = PointSetMetricType::New();
  metric->SetMovingTransform( transform );
  metric->SetFixedPointSet( fixedPoints );
  metric->SetMovingPointSet( movingPoints );
  metric->SetVirtualDomainFromImage( fixedImage );
  metric->SetCoupling( coupling );
  metric->Initialize();
  runRegistration<PointSetMetricType>(fixedPoints, movingPoints, metric, fixedImage, "ot-unbalanced.csv");

  fixedTransformedPoints =  metric->GetFixedTransformedPointSet();
  }
  clock.Stop();
  std::cout << "OT2 Metric Time : " << clock.GetTotal() << std::endl;
  }

  return EXIT_SUCCESS;
}
