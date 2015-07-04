# suppose you have created a ui file called configwin.ui
cmake_minimum_required (VERSION 2.6)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()
 
project(RegistrationApp)
 
find_package(Qt5Widgets)
 
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
 
QT5_WRAP_CPP(tst_hdr_moc configwin.h)
QT5_WRAP_UI(tst_form_hdr configwin.ui)
 
add_library(configwin ${tst_hdr_moc} ${tst_form_hdr})
qt5_use_modules(configwin Core Widgets)


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

add_library(ItkFunctionality LoadSeriesSafeAsVolume.cxx LoadSeriesSafeAsVolume.h RegisterTest2D.cxx RegisterTest2D.h Toolbox.cxx Toolbox.h ImageShort2D.h ImageShort2D.cxx itkImageToVTKImageFilter.h DeformableRegistrationTest.cpp )
add_library(VtkVisualization Visualizer.h HeatmapVisualizer.h)

add_executable(RegistrationApp Main.cpp configwin )
qt5_use_modules(RegistrationApp Core Gui Widgets)


target_link_libraries(RegistrationApp ${Glue} ${ITK_LIBRARIES} ItkFunctionality VtkVisualization ${VTK_LIBRARIES})
