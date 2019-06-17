itk_module_test()

set(OptimalTransportTests
  itkMinimalStandardRandomVariateGeneratorTest.cxx
  itkOptimalTransportImageToImageFilterTest.cxx
  )

CreateTestDriver(OptimalTransport "${OptimalTransport-Test_LIBRARIES}" "${OptimalTransportTests}")

itk_add_test(NAME itkMinimalStandardRandomVariateGeneratorTest
  COMMAND OptimalTransportTestDriver itkMinimalStandardRandomVariateGeneratorTest
  )

itk_add_test(NAME itkOptimalTransportImageToImageFilterTest
  COMMAND OptimalTransportTestDriver
    --compare
    ${ITK_TEST_OUTPUT_DIR}/itkOptimalTransportImageToImageFilterTestOutput.mha
    DATA{Baseline/itkOptimalTransportImageToImageFilterTestOutput.mha}
  itkOptimalTransportImageToImageFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkOptimalTransportImageToImageFilterTestOutput.mha
  )
