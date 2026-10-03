# fbs_genetics 0.1.0: API and host contract

`include/fbs/genetics.h` is the versioned C99 public interface; the only
implementation is `src/genetics.c`, target `fbs_genetics`, alias
`fbs::genetics`. It depends on no other FinalBuildSystems module, engine, renderer, clock,
thread, allocator, libm or global state; the library references only `memcpy`, `memmove`,
`memset` and `memcmp` (a CTest check inspects the object for `malloc`, `calloc`, `realloc` and
`free`). `FBS_GEN_VERSION` is 100 (`major*10000 + minor*100 + patch`), `FBS_GEN_SNAPSHOT_VERSION`
is 1.

The module knows no goblin words. It holds a bounded population of diploid genomes, computes trait
values from a host-authored schema, plans and commits broods, keeps a bounded pedigree with exact
kinship and inbreeding coefficients, ranks and summarizes individuals, forecasts broods and saves
portable snapshots. It never decides who breeds, never ages or kills anyone and runs on no clock of
its own: the host passes the season and maps trait values to its own rules. The goblin preset
(`examples/goblin.h`) is host data, not part of the library.

## Numbers and units

| Quantity | Type | Unit |
|---|---|---|
| Trait value, base, range, `env_sd`, `depression`, environment shift, index values | `int32_t` | milli-units: 1.000 = 1000 |
| Allele | `int8_t` | host-defined value in `[amin, amax]`, `amin > -128` |
| Dominance `h` | `uint16_t` | 0..256; 128 additive, 256 higher allele fully dominant, 0 lower allele dominant |
| Recombination fraction, mutation rate, viability floor, allele frequency, heterozygosity, `h2_q16` | `uint16_t`/`uint32_t` | Q16: 65536 = 1 |
| Kinship, inbreeding coefficient F | `uint32_t` | Q30: `FBS_GEN_Q30` = 2^30 = 1 |
| Mutation and noise scale | `uint16_t` | Q8: 256 = x1, at most `FBS_GEN_SCALE_MAX_Q8` (x16) |
| Season | `uint32_t` | host clock, monotonic |
| Serial (`fbs_gen_id`) | `uint32_t` | 1-based birth order, never reused within a context's history; 0 is never valid |

State is integer only. Doubles appear inside the library only where an output needs a variance
(`fbs_gen_stats`, forecasts, `fbs_gen_brood_preview`) or a standardized score (`fbs_gen_rank`
with `standardize = 1`, assortative pairing); they are IEEE binary64 operations built with
`-ffp-contract=off`. The tests check that identical call sequences give identical snapshot
bytes within one build; bitwise equality of floating-point summaries across compilers or
platforms (including WASM) is not established. `fbs_gen_trait_f` is a `float` convenience that never feeds
back into state.

## Model

**Loci and effects.** A schema has L loci (1..`FBS_GEN_MAX_LOCI` = 256) in map order and T traits
(1..`FBS_GEN_MAX_TRAITS` = 32). Every individual carries two alleles per locus. An effect
`(locus, trait, match, weight)` is linear, `e(a) = weight * a`, when `match == FBS_GEN_ANY`
(-128), else an indicator, `e(a) = weight if a == match`. With `hi`/`lo` the larger/smaller allele,
a locus adds `h * e(hi) + (256 - h) * e(lo)` (in 1/128 units) to its trait.

**Phenotype.** `z = clamp(base + round(sum / 128) + E - round(depression * F), lo, hi)`, rounding
half away from zero once per trait. `E` is the non-heritable part fixed at birth: the brood
environment's shift plus developmental noise. Noise is Irwin-Hall of four 16-bit uniforms scaled
to SD `env_sd * noise_scale / 256`, bounded at about ±3.46 SD. **Additive value** (`fbs_gen_additive`,
`additive_mean`) is `base + sum of e(a1) + e(a2)` over every effect (h treated as 128), no noise, no
depression or trait-range clamp (the result saturates to `int32_t`): the breeding
value when every locus is additive.

**Gametes.** A gamete walks the loci in map order: a coin picks the starting strand, and before each
later locus the strand switches with probability `recombination_q16 / 65536` (32768 = free, a new
chromosome). Each allele then mutates with probability `mutation_q16 * mutation_scale_q8 / 256 /
65536` (capped at 1) by a step uniform in 1..`mutation_step` with a random sign, clamped to
`[amin, amax]`; a clamp that leaves the allele unchanged is not counted as a mutation. A child
takes one gamete from each parent of its union: from `a` into its first strand (`a1`), from `b`
into its second (`a2`).

**Streams.** Randomness is keyed, not sequential. The child of union `u` in a brood with seed `s`
draws everything (both gametes, mutation, noise, viability) from xoshiro128\*\* 1.1 seeded by
splitmix64 from `s ^ splitmix64(u.stream)`; founder `i` of `fbs_gen_add_founders(..., seed, ...)`
uses `(seed, i)`; a pool plan draws from one stream seeded by `desc->seed`. Results therefore do
not depend on batch size, on the order unions are committed or on preview calls.

**Kinship.** `f(a, a) = (1 + F_a) / 2`; for serial `b > a` (the younger), `f(a, b) = (f(a, sire_b)
+ f(a, dam_b)) / 2`; unknown, forgotten or evicted parents contribute 0; `F(child) = f(sire, dam)`.
The window is `pedigree_depth` generations: a node that many parent steps above the query is not
expanded (its own F still counts in `f(x, x)`). Values are exact dyadic fractions in Q30 for
depths up to 8, halved with a truncating shift. The per-context memo only caches results; it never
changes a value, and it is cleared whenever a record disappears.

**Records and eviction.** A record (serial, sire, dam, born, line, F, flags, children, honours,
optionally the birth phenotype) outlives its individual so dead ancestors keep counting. When a
call needs records and the table is full, dead records born before `season - record_seasons` are
evicted lowest serial first; an evicted ancestor is afterwards unknown. If not enough records
qualify the call fails with `E_FULL` and changes nothing. `fbs_gen_forget` drops a dead record at
once.

**Lines.** A `uint32_t` per record, assigned at birth by the brood's line policy and changeable
with `fbs_gen_set_line` (the troll pen, a worker line, a house on the tablets).

## Ownership, lifetime and threading

The host owns the memory. `fbs_gen_memory_for(cfg, schema)` reports the block size (0 for an
invalid config or schema); `fbs_gen_create` lays the context out inside a block aligned for
`uint64_t` of at least that size, copying the schema (keys included), so the host's schema arrays
may be discarded afterwards. Nothing is allocated afterwards and there is no destroy: the host frees
or reuses the block. A context is not thread-safe; distinct contexts are independent. Calls that
take a non-const context but are documented as read-only (`fbs_gen_kinship`, `fbs_gen_rank`,
`fbs_gen_plan_pool`, forecasts, `fbs_gen_brood_preview`) change only caches and scratch, never
population state or snapshot bytes. Pointers returned by `fbs_gen_column` are valid until the next
call that changes the population (add, retire, forget, brood, load).

Configuration (`fbs_gen_config_default()`: 1024 living, 2048 records, depth 5, 8 seasons, 4096
memo, 2048 brood, 1024 contributors, birth traits kept):

| Field | Range | Meaning |
|---|---|---|
| `max_living` | 1..`FBS_GEN_MAX_LIVING` (2^22) | living individuals |
| `max_records` | `max_living`..`FBS_GEN_MAX_RECORDS` (2^24) | pedigree records including remembered dead |
| `pedigree_depth` | 1..8, 0 = 5 | kinship window in generations |
| `record_seasons` | 0 = 8 | dead records born within this many seasons are never evicted |
| `memo_entries` | power of two up to 2^24, 0 = 4096 | top-level kinship cache |
| `max_brood` | 1..`FBS_GEN_MAX_BROOD` (65536) | unions per plan, preview or commit |
| `max_contributors` | 1..`max_living` | contributors per pool plan |
| `record_traits` | 0 or 1 | records keep the birth phenotype (`fbs_gen_birth_traits`) |

Memory is about `max_living * (stride + 12 T + 40)` (stride = 2L rounded up to 16),
`max_records * (52 + 4T)` (4T only with `record_traits`), `max_brood * (stride + 12 T + 9)`,
`32 * max_contributors`, `16 * memo_entries` and `8 * (2^(depth+1) - 1)^2` for the per-query
kinship table (32 KB at depth 5, 2 MB at depth 8). Use `fbs_gen_memory_for` for the exact size for the complete config/schema; the
approximation omits small schema, key and alignment overheads.

## Status codes and the failure rule

Mutations are validated or staged before population changes are committed. On error,
population state and snapshot bytes are unchanged, although internal caches and scratch
may change. Outputs are unchanged with these documented exceptions: `E_TRUNCATED` size queries
report the full size through `*count`/`*bytes`; `fbs_gen_brood` and `fbs_gen_brood_preview` write
`report->failed_union` (the index of the offending union for a union error, else `n`);
`fbs_gen_brood_preview` with `cap < T` writes the first `cap` trait stats and the
complete report, then returns `E_TRUNCATED`. NULL or a readable context with an
invalid magic value returns `E_INVALID`; arbitrary invalid/dangling pointers are
not safe arguments.

| Status | Value | Meaning |
|---|---:|---|
| `FBS_GEN_OK` | 0 | success |
| `FBS_GEN_E_INVALID` | -1 | NULL required pointer, misaligned block, missing index, selfing not allowed |
| `FBS_GEN_E_RANGE` | -2 | a value outside its documented range; an invalid config at `create` |
| `FBS_GEN_E_FULL` | -3 | living or record capacity exhausted after eviction; serial space exhausted (serials never wrap); a snapshot larger than the context |
| `FBS_GEN_E_NOT_FOUND` | -4 | no record with that serial (never born, forgotten or evicted), unknown key |
| `FBS_GEN_E_DEAD` | -5 | the serial has a record but is not living |
| `FBS_GEN_E_STATE` | -6 | season backwards, forget of a living record, birth traits not kept |
| `FBS_GEN_E_SCHEMA` | -7 | bad schema; snapshot of a different structure, tuning (without retune), version or `record_traits` |
| `FBS_GEN_E_TRUNCATED` | -8 | output too small; the full size is reported |
| `FBS_GEN_E_MEMORY` | -9 | block smaller than `memory_for` reports |
| `FBS_GEN_E_CORRUPT` | -10 | snapshot bytes fail validation |

`fbs_gen_status_name` names any value (`FBS_GEN_E_UNKNOWN` otherwise).

## Calls

### Schema and context

- `fbs_gen_schema_check(schema)`: `E_INVALID` for NULL `traits`/`loci` (or `effects` with a
  non-zero count); `E_SCHEMA` for counts out of range, empty, over-long (> 32 bytes) or duplicate
  keys (traits and loci separately), `lo > base` or `base > hi`, `amin <= -128` or `amin > amax`,
  `h > 256`, `recombination_q16 > 32768` or a first locus other than 32768, a non-zero mutation
  rate with `mutation_step` outside 1..127, `founder_count` outside 1..8, a zero founder weight or
  a founder allele out of range, an effect locus or trait out of range, a `match` outside
  `[amin, amax]`; `E_RANGE` for `env_sd` outside 0..10^6, `|depression| > 2^30`, or an effect whose
  worst-case contribution (doubled per effect and added to `|base|` per trait)
  exceeds 2^30. Linear effects use `|weight| * max|allele|`; indicator effects use
  `|weight|`.
- `fbs_gen_memory_for`: 0 when the config or schema is invalid, or when the block size does not
  fit `size_t` (large configs on 32-bit and wasm32 targets; sizes are computed in 64 bits).
- `fbs_gen_create`: `E_INVALID` (NULL, block not 8-byte aligned), `E_RANGE` (config, or a size
  that does not fit `size_t`), schema statuses, `E_MEMORY`. On success the population is empty, season 0, next serial 1.
- `fbs_gen_schema_hashes`: FNV-1a 64 **structure hash** (trait and locus counts, keys in order,
  allele ranges) and **tuning hash** (trait base/lo/hi/env_sd/depression, locus h, recombination,
  mutation rate and step, founder distributions, effects in locus order). Either pointer may be NULL.
- `fbs_gen_trait_find` / `fbs_gen_locus_find`: exact byte match, `E_NOT_FOUND` otherwise.
  `fbs_gen_trait_count` / `fbs_gen_locus_count` return 0 for an invalid context.

### Season

`fbs_gen_set_season(season)`: `E_STATE` if `season` is below the current one; equal is fine.
Births record the current season; eviction eligibility uses it.

### Individuals

- `fbs_gen_add_founders(count, line, seed, out)`: `count` 0 is a no-op; `E_FULL` if living or
  records (after eviction) cannot take `count`. Founder `i` draws each allele from the locus's
  founder distribution and its noise from stream `(seed, i)`: two calls with the same seed give the
  same genomes. Flags `ALIVE | FOUNDER`, F 0, parents 0. `out` (may be NULL) receives the serials.
- `fbs_gen_add_genome(alleles, nonheritable, line, out)`: `alleles` holds `2L` bytes (`a1[L]`
  then `a2[L]`), each in its locus range (`E_RANGE`); `nonheritable` (T values, NULL for zeros)
  within ±2^29 (`E_RANGE`); `E_FULL` as above.
- `fbs_gen_retire(id)`: `E_NOT_FOUND`, `E_DEAD`. The record stays; the living slot is filled by the
  last one (slot order is deterministic and saved; outputs that order individuals tie on serial).
- `fbs_gen_forget(id)`: `E_NOT_FOUND`, `E_STATE` if living. The record is dropped; descendants
  now have an unknown parent. Clears the kinship memo.
- `fbs_gen_set_line(id, line)`: living or recorded, `E_NOT_FOUND` otherwise.
- `fbs_gen_get(id, info)`: living or recorded (`info.flags` has `FBS_GEN_ALIVE` for the living).
- `fbs_gen_traits` / `fbs_gen_additive` / `fbs_gen_nonheritable(id, out, cap)`: living only
  (`E_NOT_FOUND`, `E_DEAD`); `cap < T` is `E_TRUNCATED`.
- `fbs_gen_birth_traits`: living or recorded; `E_STATE` when `record_traits` is 0.
- `fbs_gen_genome(id, alleles, cap)`: living; `cap < 2L` is `E_TRUNCATED`.
- `fbs_gen_trait_f(id, trait)`: `z / 1000` as float, NaN on any error.
- `fbs_gen_column(trait, &ids, &values, &count)`: dense living serials and one phenotype column,
  in slot order; `E_RANGE` for a bad trait. The host's per-season sync reads these.
- `fbs_gen_records(out, cap, &count)`: every record in serial order; NULL/0 is a size query
  (`E_TRUNCATED` unless empty); `cap < count` is `E_TRUNCATED`.
- `fbs_gen_living`, `fbs_gen_record_count`, `fbs_gen_season`: counters, 0 for an invalid context.

### Kinship and pedigree

- `fbs_gen_kinship(a, b, &q30)`: both need records (`E_NOT_FOUND`); dead is fine. `a == b` gives
  `(1 + F) / 2`.
- `fbs_gen_ancestors(id, depth, out, cap, &count)`: `depth` 1..8 (`E_RANGE`, independent of
  `pedigree_depth`), `count = 2^(depth+1) - 2`, heap order (sire, dam, sire's sire, sire's dam,
  dam's sire, dam's dam, ...). An entry is the parent serial stored in the child's record, 0 for a
  founder or unknown parent; a forgotten or evicted parent keeps its serial (check it with
  `fbs_gen_get`, which then returns `E_NOT_FOUND`) and its own ancestors are 0. `E_NOT_FOUND` if `id` has no record;
  `E_TRUNCATED` when `out` is NULL or `cap < count`.

### Environment

`fbs_gen_environment_default` zeroes shifts and viability and sets both scales to 256. A brood or
forecast validates it (`E_RANGE`): scales at most 4096, `viability_count` at most 4, every one of
the 32 shifts within ±2^28, viability `trait < T` and `low <= high`. A shift is non-heritable (this
brood only). A viability rule gives survival `floor_q16` at or below `low`, 1 at or above `high`,
linear between, tested against the child's final phenotype with one draw per rule; a culled child
never enters the population and gets serial 0.

### Selection helpers

- `fbs_gen_rank(index, candidates, n, out, cap, &count)`: `index` required (`E_INVALID`), `basis`
  0..1 and `standardize` 0..1 (`E_RANGE`); `candidates` NULL or `n` 0 ranks every living
  individual; `n > max_living` is `E_RANGE`; each candidate must be living. Score: `standardize = 0`
  sums `round(weight_t * value_t / 1000)`; `standardize = 1` sums `round(weight_t * (value_t -
  mean_t) / sd_t)` over the candidate set (population SD, integer square root, 1 when 0), so a
  weight of 1000 counts one SD as 1000. Descending score, ties by ascending serial. Writes the best
  `min(cap, count)`; `*count` is the candidate count. A small `cap` is a top-k query, not an error.
- `fbs_gen_stats(trait, subset, n, out)`: over the subset (NULL/0 = everyone living): phenotype
  mean (rounded), population SD, min, max, additive mean and SD, `h2_q16 = var(A) / var(z)` capped
  at 1 (0 when var(z) is 0), mean F, `n`.
- `fbs_gen_diversity_of(subset, n, out)` and `fbs_gen_allele_frequency`: a subset longer than
  `max_living` (duplicates are allowed) is `E_RANGE`. Diversity reports loci with more than one allele, observed
  heterozygosity and expected heterozygosity `1 - sum p^2`, means over loci in Q16.
- `fbs_gen_allele_frequency(locus, allele, subset, n, &q16)`: `E_RANGE` for a bad locus.

### Breeding

**Plan.** `fbs_gen_plan_pool(desc, contributors, n, out, cap, &written)`: `E_INVALID` for NULL
arguments or an assortative plan without an index; `E_RANGE` for `children` outside
1..`max_brood`, `cap < children`, `n` outside 1..`max_contributors`, `redraws > 16`,
`kin_max_q16 > 65536`, `allow_self > 1`, `assort_k >= n` (with `n > 1`), a zero total weight;
contributor ids must be living. For each child k: draw a first parent by weight; then up to
`1 + redraws` partner draws, by weight among all contributors or, with `assort_k > 0`, by weight
among the `assort_k` contributors nearest the first parent in the assortative index (excluding
itself; ties to the lower neighbour); a disallowed self draw is skipped; stop at the first partner
with `f <= kin_max_q16` (65536 disables) and otherwise keep the least related. A child whose every
draw was a disallowed self is dropped, so `*written` can be below `children`. `union.stream =
stream_base + k` (k counts dropped children too), `union.kinship_q30` = f. The plan changes no
state.

**Commit.** `fbs_gen_brood(desc, unions, n, out_children, report)`: validation order: `desc`
(`E_INVALID`), `unions` NULL with `n > 0` (`E_INVALID`), `n > max_brood` (`E_RANGE`),
`line_policy > 2` or `allow_self > 1` (`E_RANGE`), `LINE_HIGHER` without a valid `higher_index`
(`E_INVALID`/`E_RANGE`), the environment (`E_RANGE`), then each union in order: `E_NOT_FOUND`,
`E_DEAD`, `E_INVALID` (a == b without `allow_self`), reporting its index through
`report->failed_union`; then capacity for the survivors (`E_FULL`, after eviction). The call
then evicts, appends the survivors in union order (serials ascending with the union index),
records their F, sire (`a`), dam (`b`) and line, raises each parent's `children` per surviving child
and its `honours` once per call, and writes `out_children[i]` (serial or 0 when culled) and the
report: `born`, `culled`, `mutations` (over every staged child, culled included), `evicted_records`,
`mean_f_q30` over the born, `failed_union = n`. The unions' `kinship_q30` field is ignored (F is
recomputed). Line policy: `EXPLICIT` uses `desc->line`, `FIRST` the line of `a`, `HIGHER` the line
of the parent with the higher raw (unstandardized) `higher_index` score, `a` on ties.

A plan split into chunks and committed in any order yields the same children per union (genomes,
noise, culling), provided no eviction between chunks forgets an ancestor of a later chunk's parents;
serials and honours follow the commit order.

**Preview.** `fbs_gen_brood_preview(desc, unions, n, out, cap, report)` runs the exact commit path
without committing (same validation apart from capacity) and writes per-trait stats over the
survivors in union order (`out` NULL skips the stats and `cap`); they equal `fbs_gen_stats` over the committed children in the same order,
byte for byte (tested). `report->evicted_records` is 0.

**Forecast.** `fbs_gen_forecast_union(a, b, env, out, cap, &f_q30)` (both living, `cap >= T`):
per trait it averages the four equally likely genotypes per locus, rounds the total
to milli-units, then adds the environment shift and subtracts the rounded
`depression * f(a, b)`. The variance adds per-locus
genotype variances (exact for unlinked loci; linkage changes the variance but not the mean) and
`(env_sd * noise_scale / 256)^2`; `sd` by integer square root; `p10`/`p90 = mean ∓ 1.2816 sd`.
`mean`, `p10` and `p90` are clamped to the trait range; `sd` is not. Mutation and viability are
ignored. `fbs_gen_forecast_plan(unions, n, env, out, cap)` (`n >= 1`) is the mixture: mean of
means, mean of variances plus variance of means. `tests/test_genetics_stats.c` compares a forecast
with 20,000 committed children: on an unlinked schema the mean is within 3 SE and the SD within
5%. On the linked goblin schema it prints the errors and enforces only loose bounds (mean within
5 SE plus 5 milli-units, SD within 25%): linkage and children piling up at a clamped trait range
make the forecast SD and pre-clamp mean less accurate there.

## Snapshots

`fbs_gen_snapshot_size` is exact; `fbs_gen_save(out, cap, &bytes)` reports the size and returns
`E_TRUNCATED` when `out` is NULL or `cap` is short. The same state saves to the same bytes on every
platform. Layout, little-endian throughout, no padding:

| Offset | Bytes | Field |
|---:|---:|---|
| 0 | 8 | magic `F B S G E N \0 \x01` |
| 8 | 4 | snapshot version (1) |
| 12 | 4 | API version (100) |
| 16 | 8 | structure hash |
| 24 | 8 | tuning hash |
| 32 | 4 | L |
| 36 | 4 | T |
| 40 | 4 | `record_traits` |
| 44 | 4 | season |
| 48 | 4 | next serial |
| 52 | 4 | living count |
| 56 | 4 | record count |
| 60 | 4 | reserved, 0 |
| 64 | living × (4 + 4T + 2L) | per living slot in slot order: serial, T × non-heritable part (`int32`), `a1[L]`, `a2[L]` |
| then | records × (36 + 4T if `record_traits`) | per record in serial order: serial, sire, dam, born, line, F, flags, children, honours, [T × birth phenotype] |
| end − 8 | 8 | FNV-1a 64 of every preceding byte |

Phenotypes of the living are not saved: they are recomputed from genome, non-heritable part and F
with the loading context's schema, which is what lets a save survive a balancing pass.

`fbs_gen_load(data, bytes, flags)` validates the whole image before writing anything: `E_INVALID`
(NULL), `E_RANGE` (unknown flag bits), `E_CORRUPT` (too short, magic, checksum), `E_SCHEMA`
(snapshot version, L, T or structure hash differ; tuning hash differs without
`FBS_GEN_LOAD_RETUNE`; `record_traits` differs), `E_CORRUPT` (reserved field), `E_FULL` (more
living or records than the context holds), `E_CORRUPT` for anything else: living above records,
a length mismatch, record serials zero, not strictly increasing or not below a
nonzero next serial, a sire or dam not older than its child, a birth season after the saved season, F above 1,
unknown flag bits, an alive count different from the living count, a living serial without an alive
record, a duplicate living serial, a non-heritable value beyond ±2^29, an allele out of range. Every
single-byte change is rejected (FNV-1a's per-byte step is a bijection) and leaves the context
unchanged (tested for every byte). A successful load replaces the population, season and serial
counter (0 means every serial has been used), recomputes phenotypes and clears the memo; the
API-version field at offset 12 is informational and not validated (a re-save writes the current
version); `FBS_GEN_LOAD_RETUNE` accepts a different
tuning hash with the same structure. Birth phenotypes in records are kept as saved.

## Earlier design comparison (2026-09-26)

The original lane design sketch is not included in this repository. The public
header and API contract above describe the current implementation.

Additive to the design's header sketch: `FBS_GEN_E_CORRUPT`, `FBS_GEN_MAX_RECORDS`,
`FBS_GEN_MAX_BROOD`, `FBS_GEN_MAX_VIABILITY`, `FBS_GEN_Q30`, `FBS_GEN_ENV_SD_MAX`,
`FBS_GEN_SCALE_MAX_Q8`, `fbs_gen_schema_check`, `fbs_gen_trait_count`, `fbs_gen_locus_count`,
`fbs_gen_set_line`, `fbs_gen_record_count`, `fbs_gen_nonheritable`, `fbs_gen_brood_desc.allow_self`,
`fbs_gen_brood_report.failed_union`, `fbs_gen_trait_stats.n`, the `pad` member of
`fbs_gen_scored`. Changed: `fbs_gen_rank` takes a non-const context (scratch). The snapshot
header carries `record_traits` and living entries carry no record index; records are stored sorted
and searched rather than hashed.
