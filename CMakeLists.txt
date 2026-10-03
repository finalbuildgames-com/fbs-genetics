cmake_minimum_required(VERSION 3.16)
project(FinalBuildGenetics VERSION 0.1.0 LANGUAGES C)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(_genetics_top_level OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(_genetics_top_level ON)
  include(CTest)
endif()

# Tests also need BUILD_TESTING (set by include(CTest)); a FetchContent consumer that does not
# enable testing builds the library alone. Examples default to ON only for a top-level build.
option(FBS_GENETICS_BUILD_TESTS "Build genetics C tests" ON)
option(FBS_GENETICS_BUILD_EXAMPLES "Build the genetics gathering example" ${_genetics_top_level})
option(FBS_GENETICS_BUILD_BENCHMARK "Build the genetics benchmark" OFF)
option(FBS_GENETICS_ENABLE_SANITIZERS "Enable native GCC/Clang ASan and UBSan" OFF)

add_library(fbs_genetics src/genetics.c)
add_library(fbs::genetics ALIAS fbs_genetics)
set_target_properties(fbs_genetics PROPERTIES
  EXPORT_NAME genetics OUTPUT_NAME fbs_genetics
  C_STANDARD 99 C_STANDARD_REQUIRED YES C_EXTENSIONS NO
  POSITION_INDEPENDENT_CODE ON VERSION ${PROJECT_VERSION} SOVERSION 0)
target_include_directories(fbs_genetics PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
function(fbs_genetics_c_target target)
  set_target_properties(${target} PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /fp:precise)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion
      -Wstrict-prototypes -ffp-contract=off)
  endif()
endfunction()
fbs_genetics_c_target(fbs_genetics)
if(FBS_GENETICS_ENABLE_SANITIZERS)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" AND NOT EMSCRIPTEN)
    target_compile_options(fbs_genetics PUBLIC -fsanitize=address,undefined -fno-sanitize-recover=all
      -fno-omit-frame-pointer)
    target_link_options(fbs_genetics PUBLIC -fsanitize=address,undefined)
  else()
    message(FATAL_ERROR "Genetics sanitizers require native GCC or Clang")
  endif()
endif()

set(_genetics_test_enabled OFF)
if(FBS_GENETICS_BUILD_TESTS AND BUILD_TESTING)
  set(_genetics_test_enabled ON)
  enable_testing()
  foreach(test_name genetics genetics_stats)
    add_executable(fbs_genetics_test_${test_name} tests/test_${test_name}.c)
    target_link_libraries(fbs_genetics_test_${test_name} PRIVATE fbs::genetics)
    target_include_directories(fbs_genetics_test_${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples)
    fbs_genetics_c_target(fbs_genetics_test_${test_name})
    if(NOT MSVC)
      target_link_libraries(fbs_genetics_test_${test_name} PRIVATE m)
    endif()
    add_test(NAME ${test_name} COMMAND fbs_genetics_test_${test_name})
  endforeach()
  # The library must not reference the C allocator: inspect the object's undefined symbols.
  find_program(FBS_GENETICS_NM NAMES ${CMAKE_NM} nm llvm-nm)
  if(FBS_GENETICS_NM AND NOT MSVC AND NOT EMSCRIPTEN)
    add_test(NAME genetics_no_allocation
      COMMAND ${CMAKE_COMMAND} -DNM=${FBS_GENETICS_NM} -DLIB=$<TARGET_FILE:fbs_genetics>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_no_alloc.cmake)
  endif()
endif()

if(FBS_GENETICS_BUILD_EXAMPLES)
  add_executable(fbs-genetics-gathering examples/gathering.c)
  target_link_libraries(fbs-genetics-gathering PRIVATE fbs::genetics)
  fbs_genetics_c_target(fbs-genetics-gathering)
  if(_genetics_test_enabled)
    add_test(NAME genetics_gathering COMMAND fbs-genetics-gathering 1)
  endif()
endif()

if(FBS_GENETICS_BUILD_BENCHMARK)
  add_executable(fbs-genetics-benchmark tests/benchmark.c)
  target_link_libraries(fbs-genetics-benchmark PRIVATE fbs::genetics)
  target_include_directories(fbs-genetics-benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples)
  fbs_genetics_c_target(fbs-genetics-benchmark)
endif()

install(TARGETS fbs_genetics EXPORT FinalBuildGeneticsTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES include/fbs/genetics.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fbs)
install(EXPORT FinalBuildGeneticsTargets NAMESPACE fbs::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FinalBuildGenetics)
configure_package_config_file(cmake/FinalBuildGeneticsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FinalBuildGeneticsConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FinalBuildGenetics)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/FinalBuildGeneticsConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FinalBuildGeneticsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FinalBuildGeneticsConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FinalBuildGenetics)
install(FILES examples/goblin.h DESTINATION ${CMAKE_INSTALL_DATADIR}/FinalBuildGenetics/examples)
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DATADIR}/FinalBuildGenetics/docs)

install(FILES LICENSE THIRD_PARTY_NOTICES DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/fbs-genetics)
