file(GLOB HEADERS CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/*.hpp")


add_executable(PoiseuilleFlow main.cpp)

target_link_libraries(PoiseuilleFlow PRIVATE Eigen3::Eigen)

target_include_directories(PoiseuilleFlow PUBLIC ../include)

# link python and numpy
target_link_libraries(PoiseuilleFlow
    PRIVATE
        ${PYTHON_LIBRARIES}
        Python3::NumPy
)

target_compile_features(PoiseuilleFlow PRIVATE cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )
