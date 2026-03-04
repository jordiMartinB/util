file(GLOB_RECURSE util_SRC *.cpp)
# Build object library so sources compile to object files that can be linked later.
add_library(util_objects OBJECT ${util_SRC})
set_target_properties(util_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Make headers available to other object targets that link util_objects
target_include_directories(util_objects PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
)

if(CMAKE_TESTING_ENABLED)
    add_subdirectory(tests)
endif()
