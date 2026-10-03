/* fbs/genetics.h — FinalBuildSystems heritable traits for game populations.
 *
 * C99, engine independent, no allocation, no floating point in state, no libm,
 * no global state, no threads. The host owns the memory: memory_for reports a
 * block size, create lays the context out inside it, and the block is simply
 * abandoned when the host is done (there is no destroy).
 *
 * The module holds a bounded population of diploid genomes, computes trait
 * values from a host-authored schema of loci, effects and trait ranges, plans
 * and commits broods (pairing, recombination, mutation, viability culling),
 * keeps a bounded pedigree with exact kinship and inbreeding coefficients,
 * ranks candidates, forecasts a union's brood and saves portable snapshots.
 * It never decides who breeds, never ages or kills anyone and knows nothing of
 * gameplay: the host passes the season and maps traits to its own rules.
 *
 * Numbers: trait values are int32 milli-units (1.000 = 1000); alleles are int8;
 * dominance h is 0..256 (128 = additive); probabilities and recombination
 * fractions are Q16 (65536 = 1); kinship and F are Q30 (2^30 = 1).
 *
 * Contract: modules/genetics/docs/API.md. API 0.1.0.
 */
#ifndef FBS_GENETICS_H
#define FBS_GENETICS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FBS_GEN_VERSION          100u     /* major*10000 + minor*100 + patch */
#define FBS_GEN_SNAPSHOT_VERSION 1u
#define FBS_GEN_MAX_LOCI         256u
#define FBS_GEN_MAX_TRAITS       32u
#define FBS_GEN_MAX_EFFECTS      4096u
#define FBS_GEN_MAX_LIVING       4194304u
#define FBS_GEN_MAX_RECORDS      16777216u
#define FBS_GEN_MAX_DEPTH        8u
#define FBS_GEN_MAX_BROOD        65536u
#define FBS_GEN_MAX_VIABILITY    4u
#define FBS_GEN_ANY              (-128)   /* effect.match: linear in the allele value */
#define FBS_GEN_KEY_MAX          32u
#define FBS_GEN_Q30              1073741824u
#define FBS_GEN_ENV_SD_MAX       1000000  /* milli-units */
#define FBS_GEN_SCALE_MAX_Q8     4096u    /* x16 */

/* Record flags. */
#define FBS_GEN_ALIVE   1u
#define FBS_GEN_FOUNDER 2u

typedef uint32_t fbs_gen_id;              /* serial; 0 is never valid */

typedef enum fbs_gen_status {
  FBS_GEN_OK          = 0,
  FBS_GEN_E_INVALID   = -1,  /* NULL pointer, misaligned block, bad enum, selfing not allowed */
  FBS_GEN_E_RANGE     = -2,  /* a value outside its documented range */
  FBS_GEN_E_FULL      = -3,  /* living or record capacity exhausted (after eviction) */
  FBS_GEN_E_NOT_FOUND = -4,  /* no record with that serial (never born, forgotten or evicted) */
  FBS_GEN_E_DEAD      = -5,  /* the serial has a record but is not living */
  FBS_GEN_E_STATE     = -6,  /* season went backwards, forget of a living record, ... */
  FBS_GEN_E_SCHEMA    = -7,  /* bad schema, or a snapshot of a different schema */
  FBS_GEN_E_TRUNCATED = -8,  /* output buffer too small; *count reports the full size */
  FBS_GEN_E_MEMORY    = -9,  /* block smaller than memory_for reports */
  FBS_GEN_E_CORRUPT   = -10  /* snapshot bytes fail validation */
} fbs_gen_status;

const char *fbs_gen_status_name(int status);
unsigned    fbs_gen_version(void);

/* ---- schema: authored once, copied into the block, immutable for the context's life ---- */

typedef struct fbs_gen_trait_desc {
  const char *key; size_t key_len;  /* "size", "strength", ... unique, 1..32 bytes */
  int32_t base, lo, hi;             /* milli-units; lo <= base <= hi */
  int32_t env_sd;                   /* milli-units, 0..FBS_GEN_ENV_SD_MAX; SD of developmental noise */
  int32_t depression;               /* milli-units subtracted per F = 1, may be negative */
} fbs_gen_trait_desc;

typedef struct fbs_gen_founder_allele { int8_t allele; uint16_t weight; } fbs_gen_founder_allele;

typedef struct fbs_gen_locus_desc {
  const char *key; size_t key_len;  /* "GIANT1", "SZ3" ... unique, 1..32 bytes */
  int8_t   amin, amax;              /* allele range, amin <= amax, both > -128 */
  uint16_t h;                       /* 0..256, 128 = additive, 256 = higher allele dominant */
  uint16_t recombination_q16;       /* to the previous locus, 0..32768; 32768 starts a chromosome
                                       (the first locus must be 32768) */
  uint16_t mutation_q16;            /* per allele per gamete */
  uint8_t  mutation_step;           /* max |step|, 1..127 (ignored when mutation_q16 is 0) */
  uint8_t  founder_count;           /* 1..8 entries in founder[] */
  fbs_gen_founder_allele founder[8];/* founder allele distribution, weights > 0 */
} fbs_gen_locus_desc;

/* match == FBS_GEN_ANY: e(a) = weight * a. Otherwise e(a) = weight if a == match else 0.
 * A locus contributes g = (h*e(hi) + (256-h)*e(lo)) / 128 with hi/lo the larger/smaller allele. */
typedef struct fbs_gen_effect { uint16_t locus; uint8_t trait; int8_t match; int32_t weight; } fbs_gen_effect;

typedef struct fbs_gen_schema {
  const fbs_gen_trait_desc *traits;  uint32_t trait_count;   /* 1..FBS_GEN_MAX_TRAITS */
  const fbs_gen_locus_desc *loci;    uint32_t locus_count;   /* 1..FBS_GEN_MAX_LOCI, map order */
  const fbs_gen_effect     *effects; uint32_t effect_count;  /* 0..FBS_GEN_MAX_EFFECTS */
} fbs_gen_schema;

/* ---- context ---- */

typedef struct fbs_gen_config {
  uint32_t max_living;       /* 1..FBS_GEN_MAX_LIVING */
  uint32_t max_records;      /* max_living..FBS_GEN_MAX_RECORDS; pedigree records incl. remembered dead */
  uint32_t pedigree_depth;   /* 1..8 generations of kinship recursion, 0 selects 5 */
  uint32_t record_seasons;   /* dead records born within this many seasons are never evicted, 0 selects 8 */
  uint32_t memo_entries;     /* kinship memo, power of two up to 2^24, 0 selects 4096 */
  uint32_t max_brood;        /* unions per plan/preview/commit, 1..FBS_GEN_MAX_BROOD */
  uint32_t max_contributors; /* per pool plan, 1..max_living */
  uint32_t record_traits;    /* 1: records keep the birth phenotype for ledgers */
} fbs_gen_config;
fbs_gen_config fbs_gen_config_default(void);   /* 1024 living, 2048 records, depth 5, 8 seasons,
                                                  4096 memo, 2048 brood, 1024 contributors, traits kept */

typedef struct fbs_gen fbs_gen;

/* 0 = invalid config or schema. The block must be aligned for uint64_t. */
size_t         fbs_gen_memory_for(const fbs_gen_config *cfg, const fbs_gen_schema *schema);
fbs_gen_status fbs_gen_create(const fbs_gen_config *cfg, const fbs_gen_schema *schema,
                              void *memory, size_t bytes, fbs_gen **out);
/* Detailed schema validation: FBS_GEN_OK or the reason memory_for returned 0. */
fbs_gen_status fbs_gen_schema_check(const fbs_gen_schema *schema);
fbs_gen_status fbs_gen_schema_hashes(const fbs_gen *ctx, uint64_t *structure, uint64_t *tuning);
fbs_gen_status fbs_gen_trait_find(const fbs_gen *ctx, const char *key, size_t len, uint32_t *trait);
fbs_gen_status fbs_gen_locus_find(const fbs_gen *ctx, const char *key, size_t len, uint32_t *locus);
uint32_t       fbs_gen_trait_count(const fbs_gen *ctx);
uint32_t       fbs_gen_locus_count(const fbs_gen *ctx);

/* ---- season clock: the host's ---- */
fbs_gen_status fbs_gen_set_season(fbs_gen *ctx, uint32_t season);   /* monotonic, else E_STATE */
uint32_t       fbs_gen_season(const fbs_gen *ctx);

/* ---- individuals ---- */

/* count founders, founder i drawn from each locus's founder distribution with stream (seed, i),
 * plus developmental noise. out (may be NULL) receives count serials. Atomic. */
fbs_gen_status fbs_gen_add_founders(fbs_gen *ctx, uint32_t count, uint32_t line, uint64_t seed,
                                    fbs_gen_id *out);
/* A founder with explicit alleles (2L bytes: a1[L] then a2[L]) and non-heritable part (T values,
 * may be NULL for zeros). */
fbs_gen_status fbs_gen_add_genome(fbs_gen *ctx, const int8_t *alleles, const int32_t *nonheritable,
                                  uint32_t line, fbs_gen_id *out);
fbs_gen_status fbs_gen_retire(fbs_gen *ctx, fbs_gen_id id);          /* died; the record stays */
fbs_gen_status fbs_gen_forget(fbs_gen *ctx, fbs_gen_id id);          /* drop a dead record now */
fbs_gen_status fbs_gen_set_line(fbs_gen *ctx, fbs_gen_id id, uint32_t line); /* living or recorded */
uint32_t       fbs_gen_living(const fbs_gen *ctx);
uint32_t       fbs_gen_record_count(const fbs_gen *ctx);

typedef struct fbs_gen_info {
  fbs_gen_id id, sire, dam;         /* 0 for unknown or founder */
  uint32_t born, line, children, honours, flags;   /* FBS_GEN_ALIVE, FBS_GEN_FOUNDER */
  uint32_t f_q30;                   /* own inbreeding coefficient */
} fbs_gen_info;
fbs_gen_status fbs_gen_get(const fbs_gen *ctx, fbs_gen_id id, fbs_gen_info *out);    /* living or recorded */
fbs_gen_status fbs_gen_traits(const fbs_gen *ctx, fbs_gen_id id, int32_t *z, uint32_t cap);   /* living */
/* Additive value: base + sum over loci of e(a1) + e(a2) (h = 128), no noise, no depression, no clamp. */
fbs_gen_status fbs_gen_additive(const fbs_gen *ctx, fbs_gen_id id, int32_t *a, uint32_t cap); /* living */
fbs_gen_status fbs_gen_birth_traits(const fbs_gen *ctx, fbs_gen_id id, int32_t *z, uint32_t cap); /* record */
fbs_gen_status fbs_gen_genome(const fbs_gen *ctx, fbs_gen_id id, int8_t *alleles, uint32_t cap); /* 2L, living */
fbs_gen_status fbs_gen_nonheritable(const fbs_gen *ctx, fbs_gen_id id, int32_t *e, uint32_t cap); /* living */
float          fbs_gen_trait_f(const fbs_gen *ctx, fbs_gen_id id, uint32_t trait); /* NaN on error */
/* Bulk read: dense living ids and one trait column, valid until the next mutating call. */
fbs_gen_status fbs_gen_column(const fbs_gen *ctx, uint32_t trait, const fbs_gen_id **ids,
                              const int32_t **values, uint32_t *count);
/* All records in serial order. out NULL / cap 0 is a size query. */
fbs_gen_status fbs_gen_records(const fbs_gen *ctx, fbs_gen_info *out, size_t cap, size_t *count);

/* ---- kinship ---- */
fbs_gen_status fbs_gen_kinship(fbs_gen *ctx, fbs_gen_id a, fbs_gen_id b, uint32_t *q30);
/* Heap order: sire, dam, sire's sire, sire's dam, dam's sire, ...; 0 = unknown.
 * count = 2^(depth+1) - 2, depth 1..8. */
fbs_gen_status fbs_gen_ancestors(const fbs_gen *ctx, fbs_gen_id id, uint32_t depth,
                                 fbs_gen_id *out, size_t cap, size_t *count);

/* ---- environment: brood conditions ---- */
typedef struct fbs_gen_viability {
  uint32_t trait;
  int32_t  low, high;               /* survival rises linearly from floor at <= low to 1 at >= high */
  uint16_t floor_q16;
} fbs_gen_viability;
typedef struct fbs_gen_environment {
  int32_t  shift[FBS_GEN_MAX_TRAITS];  /* non-heritable, milli-units, this brood only, |x| <= 2^28 */
  uint16_t mutation_scale_q8;          /* 256 = x1, up to FBS_GEN_SCALE_MAX_Q8 */
  uint16_t noise_scale_q8;             /* 256 = x1, up to FBS_GEN_SCALE_MAX_Q8 */
  uint32_t viability_count;            /* 0..4 */
  fbs_gen_viability viability[FBS_GEN_MAX_VIABILITY];
} fbs_gen_environment;
void fbs_gen_environment_default(fbs_gen_environment *env);

/* ---- selection helpers ---- */
typedef enum fbs_gen_basis { FBS_GEN_PHENOTYPE = 0, FBS_GEN_ADDITIVE = 1 } fbs_gen_basis;
typedef struct fbs_gen_index {
  int32_t  weight[FBS_GEN_MAX_TRAITS]; /* per trait */
  uint32_t basis;                      /* fbs_gen_basis */
  uint32_t standardize;                /* 1: score = sum weight * z-score over the candidate set;
                                          0: score = sum weight * value / 1000 */
} fbs_gen_index;
typedef struct fbs_gen_scored { fbs_gen_id id; uint32_t pad; int64_t score; } fbs_gen_scored;
/* candidates NULL / n 0 ranks every living individual. Descending score, ties by ascending serial.
 * Writes the best min(cap, n) entries; *count = n. A small cap is a top-k query, not truncation. */
fbs_gen_status fbs_gen_rank(fbs_gen *ctx, const fbs_gen_index *index,
                            const fbs_gen_id *candidates, uint32_t n,
                            fbs_gen_scored *out, size_t cap, size_t *count);

typedef struct fbs_gen_trait_stats {
  int32_t  mean, sd, min, max;         /* phenotype */
  int32_t  additive_mean, additive_sd; /* additive value */
  uint32_t h2_q16;                     /* var(additive) / var(phenotype), realized, capped at 1 */
  uint32_t mean_f_q30;
  uint32_t n;
} fbs_gen_trait_stats;
/* subset NULL / n 0 = every living individual. */
fbs_gen_status fbs_gen_stats(const fbs_gen *ctx, uint32_t trait, const fbs_gen_id *subset,
                             uint32_t n, fbs_gen_trait_stats *out);
typedef struct fbs_gen_diversity {
  uint32_t loci_polymorphic;           /* loci with more than one allele in the subset */
  uint32_t heterozygosity_q16;         /* observed, mean over loci */
  uint32_t expected_heterozygosity_q16;/* 1 - sum p_i^2, mean over loci */
} fbs_gen_diversity;
fbs_gen_status fbs_gen_diversity_of(const fbs_gen *ctx, const fbs_gen_id *subset, uint32_t n,
                                    fbs_gen_diversity *out);
fbs_gen_status fbs_gen_allele_frequency(const fbs_gen *ctx, uint32_t locus, int8_t allele,
                                        const fbs_gen_id *subset, uint32_t n, uint32_t *q16);

/* ---- breeding ---- */
typedef struct fbs_gen_union { fbs_gen_id a, b; uint32_t stream; uint32_t kinship_q30; } fbs_gen_union;
typedef struct fbs_gen_contributor { fbs_gen_id id; uint32_t weight; } fbs_gen_contributor;
typedef struct fbs_gen_pool_desc {
  uint32_t children;                   /* unions to plan, 1..max_brood */
  uint64_t seed;
  uint32_t kin_max_q16;                /* redraw while f(a,b) > this; 65536 disables */
  uint32_t redraws;                    /* 0..16 */
  uint32_t allow_self;
  uint32_t assort_k;                   /* 0 = random partner; else among the k nearest in assort_index */
  const fbs_gen_index *assort_index;
  uint32_t stream_base;                /* union.stream = stream_base + k */
} fbs_gen_pool_desc;
/* Broadcast spawning: child k draws its first parent by weight, then a partner by weight (or among
 * the assort_k nearest), redrawing up to `redraws` times while the kinship exceeds kin_max_q16 and
 * keeping the least related draw. A child whose every draw was a disallowed self is skipped, so
 * *written may be below desc->children. cap must be >= desc->children. */
fbs_gen_status fbs_gen_plan_pool(fbs_gen *ctx, const fbs_gen_pool_desc *desc,
                                 const fbs_gen_contributor *contributors, uint32_t n,
                                 fbs_gen_union *out, uint32_t cap, uint32_t *written);

typedef enum fbs_gen_line_policy {
  FBS_GEN_LINE_EXPLICIT = 0,           /* desc->line */
  FBS_GEN_LINE_FIRST    = 1,           /* the union's first parent's line */
  FBS_GEN_LINE_HIGHER   = 2            /* the parent with the higher raw higher_index score */
} fbs_gen_line_policy;
typedef struct fbs_gen_brood_desc {
  uint64_t seed;
  const fbs_gen_environment *env;      /* NULL = neutral */
  uint32_t line_policy, line;
  const fbs_gen_index *higher_index;   /* required by FBS_GEN_LINE_HIGHER */
  uint32_t allow_self;
} fbs_gen_brood_desc;
typedef struct fbs_gen_brood_report {
  uint32_t born, culled, mutations, evicted_records;
  uint32_t mean_f_q30;                 /* over the born */
  uint32_t failed_union;               /* index of the offending union on a union error, else n */
} fbs_gen_brood_report;
/* Atomic. out_children[i] (may be NULL) is the serial of union i's child, 0 when culled. Children
 * receive serials in union order. Each parent's honours count once per brood; children counts rise
 * per surviving child. The child of union u draws everything from its own stream keyed by
 * (desc->seed, u.stream), so results do not depend on batch size, order or preview calls. */
fbs_gen_status fbs_gen_brood(fbs_gen *ctx, const fbs_gen_brood_desc *desc,
                             const fbs_gen_union *unions, uint32_t n,
                             fbs_gen_id *out_children, fbs_gen_brood_report *report);

typedef struct fbs_gen_forecast { int32_t mean, sd, p10, p90; } fbs_gen_forecast;
/* Analytic brood distribution of one union, one entry per trait (cap >= trait count).
 * The analytic mean is rounded to milli-units before clamping. Per-locus variance is exact
 * for unlinked loci before rounding/clamping; configured noise variance is added. Mutation
 * and viability are ignored; p10/p90 = mean -/+ 1.2816 sd; mean, p10 and p90 are clamped
 * to the trait range. */
fbs_gen_status fbs_gen_forecast_union(fbs_gen *ctx, fbs_gen_id a, fbs_gen_id b,
                                      const fbs_gen_environment *env,
                                      fbs_gen_forecast *out, uint32_t cap, uint32_t *f_q30);
/* Mixture of the unions: mean of means, mean of variances plus variance of means. */
fbs_gen_status fbs_gen_forecast_plan(fbs_gen *ctx, const fbs_gen_union *unions, uint32_t n,
                                     const fbs_gen_environment *env,
                                     fbs_gen_forecast *out, uint32_t cap);
/* The exact commit path, nothing committed: per-trait stats over the survivors and the report
 * (evicted_records reports 0). */
fbs_gen_status fbs_gen_brood_preview(fbs_gen *ctx, const fbs_gen_brood_desc *desc,
                                     const fbs_gen_union *unions, uint32_t n,
                                     fbs_gen_trait_stats *out, uint32_t cap,
                                     fbs_gen_brood_report *report);

/* ---- snapshots ---- */
#define FBS_GEN_LOAD_RETUNE 1u   /* accept a matching structure hash with a different tuning hash */
size_t         fbs_gen_snapshot_size(const fbs_gen *ctx);
fbs_gen_status fbs_gen_save(const fbs_gen *ctx, void *out, size_t cap, size_t *bytes);
/* Validates the whole image before changing anything. Phenotypes are recomputed from genome,
 * non-heritable part and F with the context's (possibly retuned) schema. */
fbs_gen_status fbs_gen_load(fbs_gen *ctx, const void *data, size_t bytes, uint32_t flags);

#ifdef __cplusplus
}
#endif
#endif /* FBS_GENETICS_H */
