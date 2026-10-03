# Standalone genetics build

Requires CMake 3.16+ and a C99 compiler. From this repository's root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --no-tests=error
./build/fbs-genetics-gathering 1
cmake --install build --prefix ./install
```

Options: `FBS_GENETICS_BUILD_TESTS`, `FBS_GENETICS_BUILD_EXAMPLES`,
`FBS_GENETICS_BUILD_BENCHMARK` and `FBS_GENETICS_ENABLE_SANITIZERS`.
Tests also require `BUILD_TESTING`; examples default on only at top level.
The runtime has no allocation functions. Standard math is used only by tests.

An installed consumer calls `find_package(FinalBuildGenetics CONFIG REQUIRED)`
and links `fbs::genetics`, with `CMAKE_PREFIX_PATH` set to the install prefix.
FetchContent users select `https://github.com/finalbuildgames-com/fbs-genetics.git`
with a reviewed full commit as `GIT_TAG`; no `SOURCE_SUBDIR` is needed.
See the root README for a complete FetchContent example. While the repository
is private, the fetching host must have GitHub access.

The optional header-only goblin preset is in `examples/goblin.h`; a consumer
can add `${fbs_genetics_SOURCE_DIR}/examples` to its private include paths.
