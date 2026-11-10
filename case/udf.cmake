# Make headers next to channel.udf (logLawWallModel.hpp) visible to the UDF build.
target_include_directories(udf PRIVATE ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OKL PRIVATE ${CMAKE_CURRENT_LIST_DIR})
