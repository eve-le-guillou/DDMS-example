cmake_minimum_required(VERSION 3.21)

project(ttkContinuousScatterPlotCmd)

if(TARGET ttkContinuousScatterPlot)
  add_executable(${PROJECT_NAME} main.cpp)
  target_link_libraries(${PROJECT_NAME}
    PRIVATE
      ttkContinuousScatterPlot
      VTK::IOXML
    )
  set_target_properties(${PROJECT_NAME}
    PROPERTIES
      INSTALL_RPATH
        "${CMAKE_INSTALL_RPATH}"
    )
  install(
    TARGETS
      ${PROJECT_NAME}
    RUNTIME DESTINATION
      ${TTK_INSTALL_BINARY_DIR}
    )
endif()
