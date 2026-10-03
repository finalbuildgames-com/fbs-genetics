# fbs-genetics

Seeded breeding and pedigrees for game populations, in C99 with no heap allocation.

## What it does

You describe your creatures' genetics once as a schema: traits (a base value, a clamp range, developmental noise and inbreeding depression), loci (allele range, dominance, recombination with the previous locus, mutation rate and the founder allele distribution) and effects that tie loci to traits. The library then keeps a bounded population of diploid genomes and a pedigree for it. It never decides who breeds and never ages or kills anyone; your game passes the season and maps trait values to its own rules.

- `fbs_gen_memory_for` reports a block size and `fbs_gen_create` lays the context out inside a block you own.
- `fbs_gen_add_founders` draws seeded founders from the schema; `fbs_gen_add_genome` adds one with explicit alleles.
- `fbs_gen_plan_pool` pairs contributors by weight, with redraws while kinship exceeds a limit and optional assortative pairing. `fbs_gen_brood` commits the children (recombination, mutation, noise, viability culling) and returns their serials and a report. `fbs_gen_brood_preview` runs the same path without committing, and `fbs_gen_forecast_union` / `fbs_gen_forecast_plan` return an analytic mean, SD, p10 and p90 per trait.
- `fbs_gen_kinship`, `fbs_gen_get` (inbreeding coefficient F, parents, child and honour counts) and `fbs_gen_ancestors` read the pedigree. Dead individuals keep a record so they still count as ancestors.
- `fbs_gen_rank`, `fbs_gen_stats` (mean, SD, additive mean and SD, realized heritability), `fbs_gen_diversity_of`, `fbs_gen_allele_frequency` and `fbs_gen_column` cover selection and UI.
- `fbs_gen_save` / `fbs_gen_load` write and validate a portable snapshot. Loading with `FBS_GEN_LOAD_RETUNE` accepts a schema whose tuning changed but whose structure did not, and recomputes phenotypes, so a save survives a balancing pass.

Trait values are `int32_t` milli-units (1.000 = 1000), alleles are `int8_t`, probabilities are Q16 and kinship and F are Q30 (`FBS_GEN_Q30` = 1). The full contract is in [docs/API.md](docs/API.md).

## When to use it

- A game where players or the simulation breed creatures or livestock, and offspring should resemble their parents in a way that follows real inheritance (dominance, linkage, inbreeding depression).
- You need the same seed and the same calls to give the same population, and save files that are identical byte for byte.
- You want lineage UI, such as family trees and kinship warnings, or to rank breeding candidates by a weighted trait index.

## When not to use it

- Capacities are fixed at create: at most 256 loci, 32 traits, 4,194,304 living individuals, 16,777,216 pedigree records, 65,536 unions per brood and 4 viability rules per brood (`FBS_GEN_MAX_*` in the header). The defaults from `fbs_gen_config_default` are 1,024 living and 2,048 records.
- Alleles are small integers in `[amin, amax]`, not sequences. There is no notion of sex: any two living individuals can pair, the first becomes sire and the second dam.
- Kinship looks back at most `pedigree_depth` generations (1 to 8, default 5). Relatives beyond the window, and ancestors whose records were evicted or forgotten, count as unrelated.
- The forecast is exact only for unlinked loci and ignores mutation and viability. Use `fbs_gen_brood_preview` when you need the real outcome.
- A context is not thread-safe. Separate contexts are independent.
- The header has no DLL export annotations and the CMake target does not export all symbols on Windows, so link it statically there.
- This repository ships the C library only. There is no Unity, WebAssembly or engine adapter here.

## Example

Two founders, one brood of four full siblings, then each child's size and the siblings' kinship:

```c
#include <fbs/genetics.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    fbs_gen_trait_desc trait = {"size", 4, 1000, 0, 4000, 200, 0};
    fbs_gen_locus_desc locus = {"growth", 6, 0, 4, 128, 32768, 0, 0, 2, {{0, 1}, {4, 1}}};
    fbs_gen_effect effect = {0, 0, FBS_GEN_ANY, 250};
    fbs_gen_schema schema = {&trait, 1, &locus, 1, &effect, 1};
    fbs_gen_config config = fbs_gen_config_default();
    size_t bytes = fbs_gen_memory_for(&config, &schema);
    void *memory = bytes ? malloc(bytes) : NULL;
    fbs_gen *pop = NULL;
    fbs_gen_id parents[2] = {0, 0}, kids[4] = {0, 0, 0, 0};
    fbs_gen_union unions[4];
    fbs_gen_brood_desc brood;
    fbs_gen_brood_report report;
    uint32_t i, kin = 0;
    int ok;

    if (!memory) return 1;
    ok = fbs_gen_create(&config, &schema, memory, bytes, &pop) == FBS_GEN_OK &&
         fbs_gen_add_founders(pop, 2, 0, 42, parents) == FBS_GEN_OK;
    for (i = 0; i < 4; i++) {
        unions[i].a = parents[0];
        unions[i].b = parents[1];
        unions[i].stream = i;
        unions[i].kinship_q30 = 0;
    }
    memset(&brood, 0, sizeof brood); /* env NULL: neutral, nobody is culled */
    brood.seed = 7;
    ok = ok && fbs_gen_brood(pop, &brood, unions, 4, kids, &report) == FBS_GEN_OK &&
         fbs_gen_kinship(pop, kids[0], kids[1], &kin) == FBS_GEN_OK;
    if (ok) {
        for (i = 0; i < 4; i++)
            printf("child #%u size %.3f\n", (unsigned)kids[i],
                   (double)fbs_gen_trait_f(pop, kids[i], 0));
        printf("born %u, sibling kinship %.3f\n", (unsigned)report.born, (double)kin / FBS_GEN_Q30);
    }
    free(memory); /* no destroy call: the block is yours */
    return ok ? 0 : 1;
}
```

Build it with `add_executable(demo main.c)` and `target_link_libraries(demo PRIVATE fbs::genetics)`. Full siblings of unrelated founders print a kinship of 0.250.

## Build and test

Requires CMake 3.16 or newer and a C99 compiler. The library itself uses only `<string.h>`; the statistical test links libm. No third-party code is vendored.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --no-tests=error
./build/fbs-genetics-gathering 1
```

CTest runs:

- `genetics` (`tests/test_genetics.c`): config and schema limits, refused calls leaving the snapshot bytes unchanged, exact Q30 kinship (parent and offspring, full and half sibs, first cousins, selfing, Wright's full-sib series 0.25, 0.375, 0.5, 0.59375), record eviction order, allele conservation, identical runs giving identical snapshot bytes, a plan committed in scrambled chunks giving the same children, preview stats equal to commit stats, and every single-byte snapshot corruption rejected.
- `genetics_stats` (`tests/test_genetics_stats.c`, fixed seeds): 1:2:1 and 3:1 Mendelian ratios under a chi-square test, recombinant fractions within 3 SE of r, response to truncation selection, forecast mean within 3 SE and SD within 5% of 20,000 committed children on an unlinked schema, and a 12-season troll line on the goblin preset.
- `genetics_no_allocation`: `nm -u` on the library must not list `malloc`, `calloc`, `realloc`, `free`, `aligned_alloc` or `posix_memalign` (registered only when `nm` is found, not on MSVC or Emscripten).
- `genetics_gathering`: runs `fbs-genetics-gathering`, which breeds a troll line for 12 seasons on the goblin preset ([examples/goblin.h](examples/goblin.h), host data, not part of the library) and prints per-season means and a lineage tablet.

Options: `FBS_GENETICS_BUILD_TESTS` (ON, also needs `BUILD_TESTING`), `FBS_GENETICS_BUILD_EXAMPLES` (ON only for a top-level build), `FBS_GENETICS_BUILD_BENCHMARK` (OFF, builds `fbs-genetics-benchmark`, POSIX timing) and `FBS_GENETICS_ENABLE_SANITIZERS` (OFF, ASan and UBSan on GCC or Clang).

Use it from your project with FetchContent (or `add_subdirectory`) and link `fbs::genetics`:

```cmake
include(FetchContent)
FetchContent_Declare(fbs_genetics
  GIT_REPOSITORY https://github.com/finalbuildgames-com/fbs-genetics.git
  GIT_TAG <full commit hash>) # pin a reviewed commit
FetchContent_MakeAvailable(fbs_genetics)
target_link_libraries(your_game PRIVATE fbs::genetics)
```

`cmake --install build --prefix <prefix>` installs a package config, so an installed copy is found with `find_package(FinalBuildGenetics 0.1 CONFIG REQUIRED)` and `CMAKE_PREFIX_PATH` set to the prefix. See [docs/BUILD.md](docs/BUILD.md).

Engine adapters: none in this repository. Call the C API from your engine's native plugin layer.

## Design notes

- **Memory.** The host owns one block, aligned for `uint64_t`, of at least `fbs_gen_memory_for` bytes. `fbs_gen_create` copies the schema (keys included) into it, so your schema arrays can be discarded. Nothing is allocated afterwards and there is no destroy call: free or reuse the block.
- **Determinism.** Randomness is keyed, not sequential: each child draws from xoshiro128\*\* seeded by splitmix64 from the brood seed and its union's `stream`, and founder `i` from `(seed, i)`. Results therefore do not depend on batch size, commit order or preview calls. Population state is integer only. Doubles are used only for outputs (variances, standardized scores) and never feed back into state; the build uses `-ffp-contract=off`.
- **Errors.** Every call returns an `fbs_gen_status` (`fbs_gen_status_name` names it). A failed call leaves population state and snapshot bytes unchanged; documented exceptions are size queries reporting the full size and the brood report's `failed_union` index.
- **Threading.** No globals and no threads. One context must not be used from two threads at once.
- **Versioning.** `fbs_gen_version()` returns `FBS_GEN_VERSION` (100 for 0.1.0). Snapshots carry `FBS_GEN_SNAPSHOT_VERSION`, a structure hash and a tuning hash (`fbs_gen_schema_hashes`), a little-endian layout with no padding and an FNV-1a 64 checksum. `fbs_gen_load` validates the whole image before changing anything.

## License

MIT for Final Build Games' original code; see [LICENSE](LICENSE). There is no vendored third-party code; see [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES). The PRNG follows Blackman and Vigna's public-domain xoshiro128\*\* 1.1 and splitmix64.
