itk_module_test()

set(OptimalTransportTests
  itkPointSetMultiscaleOptimalTransportTest.cxx
  )

CreateTestDriver(OptimalTransport "${OptimalTransport-Test_LIBRARIES}" "${OptimalTransportTests}")

itk_add_test(NAME itkPointSetMultiscaleOptimalTransportTest
  COMMAND OptimalTransportTestDriver itkPointSetMultiscaleOptimalTransportTest

  )


