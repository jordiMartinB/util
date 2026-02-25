file(GLOB_RECURSE util_SRC *.cpp)
# Build object library so sources compile to object files that can be linked later.
add_library(util_objects OBJECT ${util_SRC})
set_target_properties(util_objects PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Make headers available to other object targets that link util_objects
target_include_directories(util_objects PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
)

# Link BZip2 if found
if (NOT BZIP2_FOUND)
    find_package(BZip2)
endif()

if (BZIP2_FOUND)
    target_compile_definitions(util_objects INTERFACE DPBUTIL_BZIP2_FOUND=${BZIP2_FOUND} )
    target_include_directories(util_objects INTERFACE ${BZIP2_INCLUDE_DIR} )
    target_link_libraries(util_objects PUBLIC ${BZIP2_LIBRARIES})
endif()

if (NOT ZLIB_FOUND)
    find_package(ZLIB)
endif()

if (ZLIB_FOUND)
    target_compile_definitions(util_objects INTERFACE DPBUTIL_ZLIB_FOUND=${ZLIB_FOUND} )
    target_include_directories(util_objects INTERFACE ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(util_objects PUBLIC ${ZLIB_LIBRARIES})
endif()

if(CMAKE_TESTING_ENABLED)
    add_subdirectory(tests)
endif()
