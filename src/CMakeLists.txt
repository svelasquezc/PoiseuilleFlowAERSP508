add_executable(PoiseuilleFlow main.cpp)

target_link_libraries(PoiseuilleFlow PRIVATE Eigen3::Eigen)

# link python and numpy
target_link_libraries(PoiseuilleFlow
    PRIVATE
        ${PYTHON_LIBRARIES}
        Python3::NumPy
)
