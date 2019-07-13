add_executable( OptimalTransportRegistrationExample OptimalTransportRegistrationExample.cxx )
 
target_link_libraries( OptimalTransportRegistrationExample ${ITK_LIBRARIES}) 

set_target_properties( OptimalTransportRegistrationExample  
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples"
)
