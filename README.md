# fbs-genetics

Heritable traits, pedigrees and deterministic population simulation.

Requires CMake 3.16+, a C99 compiler and the platform C library. No dependency
on another FinalBuildSystems module or fbs-core is needed.


```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --no-tests=error
./build/fbs-genetics-gathering 1
```

Use `add_subdirectory` or CMake FetchContent and link `fbs::genetics`:

```cmake
include(FetchContent)
FetchContent_Declare(fbs_genetics
  GIT_REPOSITORY https://github.com/finalbuildgames-com/fbs-genetics.git
  GIT_TAG main) # Pin a reviewed commit in production.
FetchContent_MakeAvailable(fbs_genetics)
target_link_libraries(your_target PRIVATE fbs::genetics)
```

Example: [gathering.c](examples/gathering.c).
Public API: [include/fbs/genetics.h](include/fbs/genetics.h).
The private FinalBuildSystems monorepo is the source of truth; releases are
curated snapshots. This repository is private pending the owner's review.

MIT for original contributions; see [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES)
for upstream licenses and attribution.
