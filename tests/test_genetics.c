/* tests/test_genetics.c — behavioural tests for fbs_genetics through the public header only.
 *
 * What it proves (docs/API.md has the contract):
 *   bounds       config and schema limits, short and misaligned blocks, capacity (E_FULL), bad
 *                arguments; every refused call leaves the snapshot bytes identical;
 *   kinship      exact Q30 values: parent-offspring, full and half sibs, first cousins, selfing,
 *                Wright's full-sib recurrence 0.25, 0.375, 0.5, 0.59375, and the documented value
 *                when the pedigree window cuts a relationship off;
 *   records      eviction order (lowest dead serial first, never inside record_seasons), forget,
 *                E_FULL when nothing is evictable, ancestors in heap order;
 *   conservation every child allele comes from the matching parent unless a mutation is counted;
 *                births grow the living count; records = living + remembered dead;
 *   determinism  identical call sequences give identical snapshot bytes; a plan committed in
 *                chunks in any order gives the same children; preview equals commit;
 *   snapshots    save/load/save identity, every single-byte corruption rejected with the state
 *                unchanged, semantic corruptions with a recomputed checksum rejected, retune only
 *                with FBS_GEN_LOAD_RETUNE, a structure change never.
 */
#include "fbs/genetics.h"
#include "goblin.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0, checks = 0;
#define CHECK(cond, ...)                                                        \
  do {                                                                          \
    checks++;                                                                   \
    if (!(cond)) {                                                              \
      failures++;                                                               \
      printf("FAIL %s:%d: ", __FILE__, __LINE__);                               \
      printf(__VA_ARGS__);                                                      \
      printf("\n");                                                             \
    }                                                                           \
  } while (0)
#define CHECK_ERR(expr, want)                                                   \
  do {                                                                          \
    int status_ = (int)(expr);                                                  \
    CHECK(status_ == (int)(want), "%s = %s, wanted %s", #expr,                  \
          fbs_gen_status_name(status_), fbs_gen_status_name((int)(want)));      \
  } while (0)
#define CHECK_OK(expr) CHECK_ERR(expr, FBS_GEN_OK)

#define Q30 1073741824u

/* ---- helpers ------------------------------------------------------------ */

static fbs_gen *make_ctx(const fbs_gen_config *cfg, const fbs_gen_schema *schema, void **mem) {
  size_t bytes = fbs_gen_memory_for(cfg, schema);
  fbs_gen *g = NULL;
  *mem = NULL;
  if (!bytes) return NULL;
  *mem = malloc(bytes);
  if (!*mem) return NULL;
  if (fbs_gen_create(cfg, schema, *mem, bytes, &g) != FBS_GEN_OK) { free(*mem); *mem = NULL; return NULL; }
  return g;
}

static uint8_t *snapshot(const fbs_gen *g, size_t *n) {
  size_t need = fbs_gen_snapshot_size(g), got = 0;
  uint8_t *p = (uint8_t *)malloc(need);
  *n = 0;
  if (p && fbs_gen_save(g, p, need, &got) != FBS_GEN_OK) { free(p); return NULL; }
  *n = got;
  return p;
}

static int same_snapshot(const fbs_gen *g, const uint8_t *img, size_t n) {
  size_t m;
  uint8_t *now = snapshot(g, &m);
  int same = now && m == n && memcmp(now, img, n) == 0;
  free(now);
  return same;
}

/* FNV-1a 64, the snapshot checksum, so tests can forge semantically bad images. */
static uint64_t fnv(const uint8_t *p, size_t n) {
  uint64_t h = 0xCBF29CE484222325ULL;
  size_t i;
  for (i = 0; i < n; i++) { h ^= p[i]; h *= 0x100000001B3ULL; }
  return h;
}
static void reseal(uint8_t *img, size_t n) {
  uint64_t h = fnv(img, n - 8u);
  unsigned i;
  for (i = 0; i < 8u; i++) img[n - 8u + i] = (uint8_t)(h >> (8u * i));
}
static void put32(uint8_t *p, uint32_t v) {
  p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8); p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
}
static uint32_t get32(const uint8_t *p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

/* One-trait schema: `loci` unlinked additive loci, alleles -8..8, founders uniform -2..2. */
typedef struct tiny_schema {
  fbs_gen_trait_desc trait;
  fbs_gen_locus_desc loci[64];
  fbs_gen_effect effects[64];
  char keys[64][12];
  fbs_gen_schema schema;
} tiny_schema;

static void tiny_build(tiny_schema *t, uint32_t loci, int32_t weight, int32_t env_sd) {
  uint32_t i, j;
  memset(t, 0, sizeof *t);
  t->trait.key = "z"; t->trait.key_len = 1;
  t->trait.base = 10000; t->trait.lo = -1000000; t->trait.hi = 1000000;
  t->trait.env_sd = env_sd;
  for (i = 0; i < loci; i++) {
    fbs_gen_locus_desc *l = &t->loci[i];
    snprintf(t->keys[i], sizeof t->keys[i], "L%u", (unsigned)i);
    l->key = t->keys[i]; l->key_len = strlen(t->keys[i]);
    l->amin = -8; l->amax = 8; l->h = 128; l->recombination_q16 = 32768;
    l->founder_count = 5;
    for (j = 0; j < 5; j++) { l->founder[j].allele = (int8_t)((int)j - 2); l->founder[j].weight = 1; }
    t->effects[i].locus = (uint16_t)i; t->effects[i].trait = 0;
    t->effects[i].match = (int8_t)FBS_GEN_ANY; t->effects[i].weight = weight;
  }
  t->schema.traits = &t->trait; t->schema.trait_count = 1;
  t->schema.loci = t->loci; t->schema.locus_count = loci;
  t->schema.effects = t->effects; t->schema.effect_count = loci;
}

static fbs_gen_config small_config(uint32_t living, uint32_t records) {
  fbs_gen_config cfg = fbs_gen_config_default();
  cfg.max_living = living; cfg.max_records = records;
  cfg.max_brood = living < 64u ? living : 64u; cfg.max_contributors = living < 64u ? living : 64u;
  return cfg;
}

static uint32_t stream_counter = 1000;

/* One child of (a, b) through a one-union brood; 0 on failure. */
static fbs_gen_id child_of(fbs_gen *g, fbs_gen_id a, fbs_gen_id b) {
  fbs_gen_brood_desc bd;
  fbs_gen_union u;
  fbs_gen_id kid = 0;
  memset(&bd, 0, sizeof bd);
  bd.seed = 77; bd.allow_self = 1;
  u.a = a; u.b = b; u.stream = stream_counter++; u.kinship_q30 = 0;
  if (fbs_gen_brood(g, &bd, &u, 1, &kid, NULL) != FBS_GEN_OK) return 0;
  return kid;
}

static fbs_gen_id founder(fbs_gen *g) {
  int8_t al[2 * 64];
  fbs_gen_id id = 0;
  memset(al, 0, sizeof al);
  if (fbs_gen_add_genome(g, al, NULL, 0, &id) != FBS_GEN_OK) return 0;
  return id;
}

static uint32_t kin(fbs_gen *g, fbs_gen_id a, fbs_gen_id b) {
  uint32_t q = 0xFFFFFFFFu;
  if (fbs_gen_kinship(g, a, b, &q) != FBS_GEN_OK) return 0xFFFFFFFFu;
  return q;
}

static uint32_t f_of(const fbs_gen *g, fbs_gen_id id) {
  fbs_gen_info info;
  if (fbs_gen_get(g, id, &info) != FBS_GEN_OK) return 0xFFFFFFFFu;
  return info.f_q30;
}

/* ---- bounds, validation and atomicity ------------------------------------ */

static void test_bounds(void) {
  tiny_schema ts;
  fbs_gen_config cfg = fbs_gen_config_default(), bad;
  fbs_gen *g = NULL;
  void *mem;
  size_t bytes, n, cnt;
  uint8_t *img;
  fbs_gen_id ids[16], kid = 0, anc[64];
  fbs_gen_union u[4];
  fbs_gen_brood_desc bd;
  fbs_gen_brood_report rep;
  fbs_gen_pool_desc pd;
  fbs_gen_contributor con[4];
  fbs_gen_environment env;
  fbs_gen_info info;
  fbs_gen_index ix;
  fbs_gen_scored sc[16];
  fbs_gen_trait_stats stats;
  fbs_gen_forecast fc[2];
  uint32_t q, written;
  int32_t z[4];
  int8_t al[8];
  const fbs_gen_id *col_ids; const int32_t *col_z;

  CHECK(fbs_gen_version() == FBS_GEN_VERSION, "version");
  CHECK(strcmp(fbs_gen_status_name(FBS_GEN_E_CORRUPT), "FBS_GEN_E_CORRUPT") == 0, "status name");
  tiny_build(&ts, 2, 100, 0);
  CHECK(fbs_gen_memory_for(&cfg, &ts.schema) > 0u, "default config is valid");
  CHECK(fbs_gen_memory_for(NULL, &ts.schema) == 0u, "NULL config");
  CHECK(fbs_gen_memory_for(&cfg, NULL) == 0u, "NULL schema");
  bad = cfg; bad.max_living = 0; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "living 0");
  bad = cfg; bad.max_living = FBS_GEN_MAX_LIVING + 1u; bad.max_records = bad.max_living;
  CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "living > max");
  bad = cfg; bad.max_records = cfg.max_living - 1u; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "records < living");
  bad = cfg; bad.max_records = FBS_GEN_MAX_RECORDS + 1u; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "records > max");
  bad = cfg; bad.pedigree_depth = 9; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "depth 9");
  bad = cfg; bad.memo_entries = 1000; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "memo not pow2");
  bad = cfg; bad.memo_entries = 1u << 25; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "memo too big");
  bad = cfg; bad.max_brood = 0; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "brood 0");
  bad = cfg; bad.max_brood = FBS_GEN_MAX_BROOD + 1u; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "brood > max");
  bad = cfg; bad.max_contributors = 0; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "contributors 0");
  bad = cfg; bad.max_contributors = cfg.max_living + 1u; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "contributors");
  bad = cfg; bad.record_traits = 2; CHECK(fbs_gen_memory_for(&bad, &ts.schema) == 0u, "record_traits 2");
  bad = cfg; bad.pedigree_depth = 0; bad.record_seasons = 0; bad.memo_entries = 0;
  CHECK(fbs_gen_memory_for(&bad, &ts.schema) == fbs_gen_memory_for(&cfg, &ts.schema), "zeros select defaults");

  /* schema errors */
  CHECK_ERR(fbs_gen_schema_check(NULL), FBS_GEN_E_INVALID);
  CHECK_OK(fbs_gen_schema_check(&ts.schema));
  ts.loci[1].key = "L0"; ts.loci[1].key_len = 2;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* duplicate locus key */
  tiny_build(&ts, 2, 100, 0); ts.loci[1].key_len = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* empty key */
  tiny_build(&ts, 2, 100, 0); ts.loci[1].key_len = FBS_GEN_KEY_MAX + 1u;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* long key */
  tiny_build(&ts, 2, 100, 0); ts.loci[0].h = 257;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.loci[0].recombination_q16 = 100;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* first locus starts a chromosome */
  tiny_build(&ts, 2, 100, 0); ts.loci[1].recombination_q16 = 32769;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.loci[1].amin = 9;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* amin > amax */
  tiny_build(&ts, 2, 100, 0); ts.loci[1].amin = -128;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* -128 is FBS_GEN_ANY */
  tiny_build(&ts, 2, 100, 0); ts.loci[1].mutation_q16 = 10; ts.loci[1].mutation_step = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.loci[1].founder_count = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.loci[1].founder[2].weight = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.loci[1].founder[2].allele = 9;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.effects[1].locus = 7;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.effects[1].trait = 1;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.effects[0].match = 9;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* match outside the allele range */
  tiny_build(&ts, 2, 100, 0); ts.trait.lo = 20000;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);           /* lo > base */
  tiny_build(&ts, 2, 100, 0); ts.trait.env_sd = -1;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_RANGE);
  tiny_build(&ts, 2, 100, 0); ts.schema.trait_count = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.schema.locus_count = 0;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_SCHEMA);
  tiny_build(&ts, 2, 100, 0); ts.schema.effects = NULL;
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_INVALID);
  tiny_build(&ts, 2, 0x7FFFFFFF / 8, 0);
  CHECK_ERR(fbs_gen_schema_check(&ts.schema), FBS_GEN_E_RANGE);            /* overflow proof */
  tiny_build(&ts, 2, 100, 0);
  CHECK_ERR(fbs_gen_create(&cfg, NULL, NULL, 0, &g), FBS_GEN_E_INVALID);

  /* short and misaligned blocks */
  cfg = small_config(8, 12); cfg.record_seasons = 2;
  bytes = fbs_gen_memory_for(&cfg, &ts.schema);
  mem = malloc(bytes + 16u);
  CHECK_ERR(fbs_gen_create(&cfg, &ts.schema, mem, bytes - 1u, &g), FBS_GEN_E_MEMORY);
  CHECK_ERR(fbs_gen_create(&cfg, &ts.schema, (uint8_t *)mem + 4, bytes, &g), FBS_GEN_E_INVALID);
  bad = cfg; bad.max_living = 0;
  CHECK_ERR(fbs_gen_create(&bad, &ts.schema, mem, bytes, &g), FBS_GEN_E_RANGE);
  CHECK_OK(fbs_gen_create(&cfg, &ts.schema, mem, bytes, &g));
  CHECK(fbs_gen_trait_count(g) == 1u && fbs_gen_locus_count(g) == 2u, "counts");
  CHECK(fbs_gen_trait_count(NULL) == 0u && fbs_gen_living(NULL) == 0u, "NULL ctx");
  CHECK_OK(fbs_gen_trait_find(g, "z", 1, &q)); CHECK(q == 0u, "trait_find");
  CHECK_ERR(fbs_gen_trait_find(g, "y", 1, &q), FBS_GEN_E_NOT_FOUND);
  CHECK_OK(fbs_gen_locus_find(g, "L1", 2, &q)); CHECK(q == 1u, "locus_find");
  CHECK_ERR(fbs_gen_locus_find(g, "L2", 2, &q), FBS_GEN_E_NOT_FOUND);

  /* capacity: 8 living */
  CHECK_OK(fbs_gen_add_founders(g, 6, 0, 1, ids));
  CHECK(fbs_gen_living(g) == 6u && fbs_gen_record_count(g) == 6u, "six founders");
  img = snapshot(g, &n);
  CHECK_ERR(fbs_gen_add_founders(g, 3, 0, 2, ids + 6), FBS_GEN_E_FULL);
  CHECK(same_snapshot(g, img, n), "refused add_founders changed state");
  memset(al, 0, sizeof al); al[0] = 9;
  CHECK_ERR(fbs_gen_add_genome(g, al, NULL, 0, &kid), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_add_genome(g, NULL, NULL, 0, &kid), FBS_GEN_E_INVALID);
  al[0] = 0; z[0] = 1 << 30;
  CHECK_ERR(fbs_gen_add_genome(g, al, z, 0, &kid), FBS_GEN_E_RANGE);
  CHECK(same_snapshot(g, img, n), "refused add_genome changed state");

  memset(&bd, 0, sizeof bd);
  bd.seed = 5;
  u[0].a = ids[0]; u[0].b = ids[1]; u[0].stream = 0; u[0].kinship_q30 = 0;
  u[1].a = ids[2]; u[1].b = ids[3]; u[1].stream = 1; u[1].kinship_q30 = 0;
  u[2].a = ids[4]; u[2].b = ids[5]; u[2].stream = 2; u[2].kinship_q30 = 0;
  memset(&rep, 0xAB, sizeof rep);
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 3, NULL, &rep), FBS_GEN_E_FULL);      /* 6 + 3 > 8 */
  CHECK(rep.failed_union == 3u, "E_FULL failed_union %u", rep.failed_union);
  CHECK(same_snapshot(g, img, n), "refused brood changed state");
  u[1].b = 999;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_NOT_FOUND);
  CHECK(rep.failed_union == 1u, "failed union index %u", rep.failed_union);
  u[1].b = ids[2];
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_INVALID);  /* selfing not allowed */
  CHECK(rep.failed_union == 1u, "selfing index %u", rep.failed_union);
  u[1].b = ids[3];
  bd.line_policy = 3;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_RANGE);
  CHECK(rep.failed_union == 2u, "non-union error reports n, got %u", rep.failed_union);
  bd.line_policy = FBS_GEN_LINE_HIGHER;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_INVALID);  /* HIGHER needs an index */
  bd.line_policy = FBS_GEN_LINE_EXPLICIT;
  fbs_gen_environment_default(&env);
  env.mutation_scale_q8 = FBS_GEN_SCALE_MAX_Q8 + 1u; bd.env = &env;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_RANGE);
  fbs_gen_environment_default(&env); env.viability_count = 5;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_RANGE);
  fbs_gen_environment_default(&env); env.viability_count = 1; env.viability[0].trait = 1;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_RANGE);
  fbs_gen_environment_default(&env); env.shift[0] = (1 << 28) + 1;
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 2, NULL, &rep), FBS_GEN_E_RANGE);
  bd.env = NULL;
  CHECK_ERR(fbs_gen_brood(g, NULL, u, 2, NULL, &rep), FBS_GEN_E_INVALID);
  CHECK_ERR(fbs_gen_brood(g, &bd, NULL, 2, NULL, &rep), FBS_GEN_E_INVALID);
  CHECK_ERR(fbs_gen_brood(g, &bd, u, cfg.max_brood + 1u, NULL, &rep), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_brood_preview(g, &bd, u, cfg.max_brood + 1u, NULL, 0, &rep), FBS_GEN_E_RANGE);
  CHECK(same_snapshot(g, img, n), "refused broods changed state");

  /* plan_pool argument checks */
  memset(&pd, 0, sizeof pd);
  pd.children = 2; pd.seed = 3; pd.kin_max_q16 = 65536;
  con[0].id = ids[0]; con[0].weight = 1; con[1].id = ids[1]; con[1].weight = 1;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 1, &written), FBS_GEN_E_RANGE);  /* cap < children */
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 0, u, 4, &written), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, cfg.max_contributors + 1u, u, 4, &written), FBS_GEN_E_RANGE);
  pd.redraws = 17;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_RANGE);
  pd.redraws = 0; pd.kin_max_q16 = 65537;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_RANGE);
  pd.kin_max_q16 = 65536; pd.children = 0;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_RANGE);
  pd.children = 2; con[1].weight = 0; con[0].weight = 0;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_RANGE); /* zero total weight */
  con[0].weight = 1; con[1].weight = 1; con[1].id = 999;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_NOT_FOUND);
  con[1].id = ids[1]; pd.assort_k = 1;
  CHECK_ERR(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written), FBS_GEN_E_INVALID); /* no assort index */
  pd.assort_k = 0; pd.redraws = 16;
  CHECK_OK(fbs_gen_plan_pool(g, &pd, con, 2, u, 4, &written));
  CHECK(written == 2u && u[0].a != u[0].b, "plan without selfing");
  /* a single contributor and no selfing: every child skipped */
  CHECK_OK(fbs_gen_plan_pool(g, &pd, con, 1, u, 4, &written));
  CHECK(written == 0u, "single contributor, no selfing: %u", written);
  pd.allow_self = 1;
  CHECK_OK(fbs_gen_plan_pool(g, &pd, con, 1, u, 4, &written));
  CHECK(written == 2u && u[0].a == u[0].b, "single contributor with selfing");
  CHECK(same_snapshot(g, img, n), "plan_pool changed state");

  /* readers */
  CHECK_ERR(fbs_gen_traits(g, ids[0], z, 0), FBS_GEN_E_TRUNCATED);
  CHECK_ERR(fbs_gen_traits(g, 999, z, 4), FBS_GEN_E_NOT_FOUND);
  CHECK_ERR(fbs_gen_genome(g, ids[0], al, 3), FBS_GEN_E_TRUNCATED);
  CHECK_OK(fbs_gen_genome(g, ids[0], al, 4));
  CHECK_ERR(fbs_gen_column(g, 1, &col_ids, &col_z, &q), FBS_GEN_E_RANGE);
  CHECK_OK(fbs_gen_column(g, 0, &col_ids, &col_z, &q));
  CHECK(q == 6u, "column count");
  CHECK_ERR(fbs_gen_records(g, NULL, 0, &cnt), FBS_GEN_E_TRUNCATED);
  CHECK(cnt == 6u, "records size query %u", (unsigned)cnt);
  CHECK(fbs_gen_trait_f(g, 999, 0) != fbs_gen_trait_f(g, 999, 0), "NaN on error");
  CHECK_ERR(fbs_gen_ancestors(g, ids[0], 0, anc, 64, &cnt), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_ancestors(g, ids[0], 3, anc, 13, &cnt), FBS_GEN_E_TRUNCATED);
  CHECK(cnt == 14u, "ancestors count %u", (unsigned)cnt);
  CHECK_ERR(fbs_gen_stats(g, 1, NULL, 0, &stats), FBS_GEN_E_RANGE);
  {
    fbs_gen_id many[9];
    fbs_gen_diversity dv;
    uint32_t k;
    for (k = 0; k < 9u; k++) many[k] = ids[0];
    CHECK_ERR(fbs_gen_diversity_of(g, many, 9, &dv), FBS_GEN_E_RANGE);          /* n > max_living */
    CHECK_ERR(fbs_gen_allele_frequency(g, 0, 0, many, 9, &q), FBS_GEN_E_RANGE);
    CHECK_OK(fbs_gen_diversity_of(g, many, 8, &dv));
    CHECK(dv.expected_heterozygosity_q16 <= 65536u, "diversity of duplicates");
  }
  CHECK_ERR(fbs_gen_forecast_union(g, ids[0], ids[1], NULL, fc, 0, NULL), FBS_GEN_E_TRUNCATED);
  CHECK_ERR(fbs_gen_forecast_union(g, ids[0], 999, NULL, fc, 2, NULL), FBS_GEN_E_NOT_FOUND);
  memset(&ix, 0, sizeof ix); ix.basis = 2;
  CHECK_ERR(fbs_gen_rank(g, &ix, NULL, 0, sc, 16, &cnt), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_rank(g, NULL, NULL, 0, sc, 16, &cnt), FBS_GEN_E_INVALID);
  ix.basis = FBS_GEN_PHENOTYPE; ix.weight[0] = 1000;
  CHECK_OK(fbs_gen_rank(g, &ix, NULL, 0, sc, 2, &cnt));
  CHECK(cnt == 6u && sc[0].score >= sc[1].score, "top-k rank");

  /* season, retire, forget */
  CHECK_OK(fbs_gen_set_season(g, 3));
  CHECK_ERR(fbs_gen_set_season(g, 2), FBS_GEN_E_STATE);
  CHECK(fbs_gen_season(g) == 3u, "season");
  CHECK_ERR(fbs_gen_forget(g, ids[0]), FBS_GEN_E_STATE);                    /* living */
  CHECK_ERR(fbs_gen_retire(g, 999), FBS_GEN_E_NOT_FOUND);
  CHECK_OK(fbs_gen_retire(g, ids[0]));
  CHECK_ERR(fbs_gen_retire(g, ids[0]), FBS_GEN_E_DEAD);
  CHECK_ERR(fbs_gen_traits(g, ids[0], z, 4), FBS_GEN_E_DEAD);
  CHECK_OK(fbs_gen_get(g, ids[0], &info));
  CHECK(!(info.flags & FBS_GEN_ALIVE) && (info.flags & FBS_GEN_FOUNDER), "retired founder flags");
  u[0].a = ids[0]; u[0].b = ids[1];
  CHECK_ERR(fbs_gen_brood(g, &bd, u, 1, NULL, &rep), FBS_GEN_E_DEAD);
  CHECK(rep.failed_union == 0u, "dead union index");
  CHECK_OK(fbs_gen_forget(g, ids[0]));
  CHECK_ERR(fbs_gen_get(g, ids[0], &info), FBS_GEN_E_NOT_FOUND);
  CHECK(fbs_gen_living(g) == 5u && fbs_gen_record_count(g) == 5u, "after forget");
  free(img);

  /* load refusals leave state unchanged */
  img = snapshot(g, &n);
  CHECK_ERR(fbs_gen_load(g, img, n, 2u), FBS_GEN_E_RANGE);
  CHECK_ERR(fbs_gen_load(g, img, n - 1u, 0), FBS_GEN_E_CORRUPT);
  CHECK_ERR(fbs_gen_load(g, img, 10, 0), FBS_GEN_E_CORRUPT);
  CHECK_ERR(fbs_gen_load(g, NULL, n, 0), FBS_GEN_E_INVALID);
  CHECK(same_snapshot(g, img, n), "refused loads changed state");
  free(img);
  free(mem);

  /* record capacity: records full of recent dead cannot be evicted */
  cfg = small_config(4, 4); cfg.record_seasons = 5;
  g = make_ctx(&cfg, &ts.schema, &mem);
  CHECK(g != NULL, "records ctx");
  if (g) {
    CHECK_OK(fbs_gen_add_founders(g, 4, 0, 9, ids));
    CHECK_OK(fbs_gen_retire(g, ids[3]));
    img = snapshot(g, &n);
    CHECK_ERR(fbs_gen_add_founders(g, 1, 0, 1, &kid), FBS_GEN_E_FULL);
    u[0].a = ids[0]; u[0].b = ids[1]; u[0].stream = 0;
    CHECK_ERR(fbs_gen_brood(g, &bd, u, 1, NULL, &rep), FBS_GEN_E_FULL);
    CHECK(same_snapshot(g, img, n), "refused record capacity changed state");
    free(img);
    free(mem);
  }
}

/* ---- kinship -------------------------------------------------------------- */

static void test_kinship(void) {
  tiny_schema ts;
  fbs_gen_config cfg;
  fbs_gen *g;
  void *mem;
  fbs_gen_id A, B, C, D, S1, S2, H, P1, P2, K1, K2, X, Y, Z;
  fbs_gen_id sib[2], anc[14];
  size_t cnt;
  uint32_t t;
  static const uint32_t wright[4] = { Q30 / 4u, Q30 / 8u * 3u, Q30 / 2u, Q30 / 32u * 19u };

  tiny_build(&ts, 4, 100, 0);
  cfg = small_config(256, 512); cfg.pedigree_depth = 8;
  g = make_ctx(&cfg, &ts.schema, &mem);
  CHECK(g != NULL, "kinship ctx");
  if (!g) return;
  X = founder(g);                                          /* forgotten later: shifts every record */
  A = founder(g); B = founder(g); C = founder(g); D = founder(g);
  CHECK(kin(g, A, A) == Q30 / 2u, "self kinship of a founder %u", kin(g, A, A));
  CHECK(kin(g, A, B) == 0u, "unrelated founders");
  S1 = child_of(g, A, B); S2 = child_of(g, A, B);          /* full sibs */
  H = child_of(g, A, C);                                   /* half sib of S1 through A */
  CHECK(S1 && S2 && H, "children");
  CHECK(kin(g, A, S1) == Q30 / 4u, "parent-offspring %u", kin(g, A, S1));
  CHECK(kin(g, S1, A) == Q30 / 4u, "parent-offspring symmetric");
  CHECK(kin(g, S1, S2) == Q30 / 4u, "full sibs %u", kin(g, S1, S2));
  CHECK(kin(g, S1, H) == Q30 / 8u, "half sibs %u", kin(g, S1, H));
  CHECK(f_of(g, S1) == 0u, "F of an outbred child");
  P1 = child_of(g, S1, C); P2 = child_of(g, S2, D);        /* first cousins */
  CHECK(kin(g, P1, P2) == Q30 / 16u, "first cousins %u", kin(g, P1, P2));
  K1 = child_of(g, P1, P2);
  CHECK(f_of(g, K1) == Q30 / 16u, "child of first cousins F %u", f_of(g, K1));
  CHECK(kin(g, K1, K1) == (Q30 + Q30 / 16u) / 2u, "self kinship includes F");
  X = child_of(g, A, A);                                   /* selfing */
  CHECK(f_of(g, X) == Q30 / 2u, "selfing F %u", f_of(g, X));
  K2 = child_of(g, S1, S2);                                 /* full-sib mating */
  CHECK(f_of(g, K2) == Q30 / 4u, "full-sib child F %u", f_of(g, K2));

  /* ancestors in heap order */
  CHECK_OK(fbs_gen_ancestors(g, K1, 3, anc, 14, &cnt));
  CHECK(cnt == 14u && anc[0] == P1 && anc[1] == P2 && anc[2] == S1 && anc[3] == C && anc[4] == S2 && anc[5] == D,
        "ancestors level 1-2");
  CHECK(anc[6] == A && anc[7] == B && anc[8] == 0u && anc[9] == 0u && anc[10] == A && anc[11] == B, "ancestors level 3");

  /* Wright: repeated full-sib mating F_t = (1 + 2F_{t-1} + F_{t-2}) / 4 */
  sib[0] = child_of(g, C, D); sib[1] = child_of(g, C, D);
  for (t = 0; t < 4u; t++) {
    Y = child_of(g, sib[0], sib[1]);
    Z = child_of(g, sib[0], sib[1]);
    CHECK(Y && Z, "wright children");
    CHECK(f_of(g, Y) == wright[t] && f_of(g, Z) == wright[t], "Wright generation %u: F %u, wanted %u",
          (unsigned)t + 1u, f_of(g, Y), wright[t]);
    sib[0] = Y; sib[1] = Z;
  }
  /* forgetting an unrelated record compacts the table; every value stays */
  CHECK_OK(fbs_gen_retire(g, X));
  CHECK_OK(fbs_gen_forget(g, X));
  CHECK(kin(g, P1, P2) == Q30 / 16u && kin(g, S1, H) == Q30 / 8u && kin(g, Y, Z) == Q30 / 64u * 43u, /* F5 = (1 + 2 F4 + F3) / 4 */
        "kinship after compaction: cousins %u half sibs %u", kin(g, P1, P2), kin(g, S1, H));
  /* forgetting an ancestor makes it unknown: full sibs through A and B keep only B's share */
  CHECK_OK(fbs_gen_retire(g, A));
  CHECK_OK(fbs_gen_forget(g, A));
  CHECK(kin(g, S1, S2) == Q30 / 8u, "full sibs with a forgotten sire %u", kin(g, S1, S2));
  CHECK(kin(g, P1, P2) == Q30 / 32u, "cousins with a forgotten great-grandparent %u", kin(g, P1, P2));
  {
    /* a snapshot round trip keeps the pedigree walk */
    size_t n;
    uint8_t *img = snapshot(g, &n);
    void *mem2;
    fbs_gen *h = make_ctx(&cfg, &ts.schema, &mem2);
    CHECK(img && h, "round trip ctx");
    if (img && h) {
      CHECK_OK(fbs_gen_load(h, img, n, 0));
      CHECK(kin(h, S1, S2) == Q30 / 8u && kin(h, S1, H) == kin(g, S1, H) && kin(h, Y, Z) == kin(g, Y, Z),
            "kinship after load");
      free(mem2);
    }
    free(img);
  }
  free(mem);

  /* the window: depth 1 sees parents only, so first cousins read as unrelated */
  cfg.pedigree_depth = 1;
  g = make_ctx(&cfg, &ts.schema, &mem);
  CHECK(g != NULL, "window ctx");
  if (!g) return;
  A = founder(g); B = founder(g); C = founder(g); D = founder(g);
  S1 = child_of(g, A, B); S2 = child_of(g, A, B);
  P1 = child_of(g, S1, C); P2 = child_of(g, S2, D);
  CHECK(kin(g, S1, S2) == Q30 / 4u, "depth 1 full sibs %u", kin(g, S1, S2));
  CHECK(kin(g, A, S1) == Q30 / 4u, "depth 1 parent-offspring");
  CHECK(kin(g, P1, P2) == 0u, "depth 1 cuts first cousins to 0, got %u", kin(g, P1, P2));
  CHECK(kin(g, P1, A) == 0u, "depth 1 cuts grandparent to 0, got %u", kin(g, P1, A));
  free(mem);
  cfg.pedigree_depth = 2;
  g = make_ctx(&cfg, &ts.schema, &mem);
  if (!g) return;
  A = founder(g); B = founder(g); C = founder(g); D = founder(g);
  S1 = child_of(g, A, B); S2 = child_of(g, A, B);
  P1 = child_of(g, S1, C); P2 = child_of(g, S2, D);
  CHECK(kin(g, P1, P2) == Q30 / 16u, "depth 2 first cousins %u", kin(g, P1, P2));
  CHECK(kin(g, P1, A) == Q30 / 8u, "depth 2 grandparent %u", kin(g, P1, A));
  free(mem);
}

/* ---- records: eviction and forget ----------------------------------------- */

static void test_eviction(void) {
  tiny_schema ts;
  fbs_gen_config cfg;
  fbs_gen *g;
  void *mem;
  fbs_gen_id ids[8], kid;
  fbs_gen_info info[8];
  fbs_gen_brood_desc bd;
  fbs_gen_brood_report rep;
  fbs_gen_union u;
  size_t cnt;

  tiny_build(&ts, 2, 100, 0);
  cfg = small_config(5, 6); cfg.record_seasons = 2;
  g = make_ctx(&cfg, &ts.schema, &mem);
  CHECK(g != NULL, "eviction ctx");
  if (!g) return;
  CHECK_OK(fbs_gen_add_founders(g, 4, 0, 3, ids));          /* season 0 */
  CHECK_OK(fbs_gen_set_season(g, 1));
  CHECK_OK(fbs_gen_add_founders(g, 0, 0, 3, NULL));
  CHECK_OK(fbs_gen_retire(g, ids[1]));
  CHECK_OK(fbs_gen_retire(g, ids[0]));
  memset(&bd, 0, sizeof bd); bd.seed = 1;
  u.a = ids[2]; u.b = ids[3]; u.stream = 0; u.kinship_q30 = 0;
  CHECK_OK(fbs_gen_brood(g, &bd, &u, 1, &kid, &rep));       /* record 5 of 6 */
  u.stream = 1;
  CHECK_OK(fbs_gen_brood(g, &bd, &u, 1, &kid, &rep));       /* record 6 of 6 */
  u.stream = 2;
  CHECK_ERR(fbs_gen_brood(g, &bd, &u, 1, &kid, &rep), FBS_GEN_E_FULL);  /* dead but too young */
  CHECK_OK(fbs_gen_set_season(g, 3));                        /* born 0 < 3 - 2: evictable */
  CHECK_OK(fbs_gen_brood(g, &bd, &u, 1, &kid, &rep));
  CHECK(rep.evicted_records == 1u, "one record evicted, got %u", rep.evicted_records);
  CHECK_ERR(fbs_gen_get(g, ids[0], info), FBS_GEN_E_NOT_FOUND);         /* lowest dead serial first */
  CHECK_OK(fbs_gen_get(g, ids[1], info));
  CHECK_OK(fbs_gen_records(g, info, 8, &cnt));
  CHECK(cnt == 6u && info[0].id == ids[1] && info[1].id == ids[2], "records in serial order");
  CHECK(fbs_gen_living(g) == 5u, "living %u", fbs_gen_living(g));
  CHECK_OK(fbs_gen_get(g, ids[2], info));
  CHECK(info[0].children == 3u && info[0].honours == 3u, "children %u honours %u", info[0].children, info[0].honours);
  {
    size_t k;
    fbs_gen_id sibs[3];
    uint32_t q, ns = 0;
    CHECK_OK(fbs_gen_records(g, info, 8, &cnt));
    for (k = 0; k < cnt; k++) if (info[k].sire == ids[2]) sibs[ns++] = info[k].id;
    CHECK(ns == 3u, "three children recorded");
    if (ns == 3u) {
      CHECK_OK(fbs_gen_kinship(g, sibs[0], sibs[2], &q));
      CHECK(q == Q30 / 4u, "full sibs across an eviction %u", q);
    }
  }
  free(mem);
}

/* ---- conservation and preview/commit on the goblin schema ---------------- */

typedef struct gob_ctx { goblin_schema_data sd; fbs_gen *g; void *mem; } gob_ctx;

static int gob_open(gob_ctx *c, uint32_t living, uint32_t brood) {
  fbs_gen_config cfg = fbs_gen_config_default();
  goblin_schema_build(&c->sd, 40);
  cfg.max_living = living; cfg.max_records = living * 2u; cfg.max_brood = brood; cfg.max_contributors = 512;
  if (cfg.max_contributors > living) cfg.max_contributors = living;
  c->g = make_ctx(&cfg, &c->sd.schema, &c->mem);
  return c->g != NULL;
}

static void test_conservation(void) {
  gob_ctx *c = (gob_ctx *)malloc(sizeof *c);
  fbs_gen_id ids[200], kids[400];
  fbs_gen_contributor con[200];
  fbs_gen_union un[400];
  fbs_gen_pool_desc pd;
  fbs_gen_brood_desc bd;
  fbs_gen_brood_report rep;
  fbs_gen_environment env;
  uint32_t i, l, written, L, mism, round, living0, records0;
  int8_t ga[96], gb[96], gk[96];
  size_t cnt;
  if (!c || !gob_open(c, 2048, 512)) { CHECK(0, "goblin ctx"); free(c); return; }
  L = fbs_gen_locus_count(c->g);
  CHECK(L == 48u, "goblin loci %u", L);
  CHECK_OK(fbs_gen_add_founders(c->g, 200, 0, 11, ids));
  for (i = 0; i < 200u; i++) { con[i].id = ids[i]; con[i].weight = 1u + i % 3u; }
  for (round = 0; round < 2u; round++) {
    memset(&pd, 0, sizeof pd);
    pd.children = 400; pd.seed = 90u + round; pd.kin_max_q16 = 65536;
    CHECK_OK(fbs_gen_plan_pool(c->g, &pd, con, 200, un, 400, &written));
    CHECK(written > 380u && written <= 400u, "planned %u (a draw of self without redraws is skipped)", written);
    memset(&bd, 0, sizeof bd);
    bd.seed = 1234u + round;
    fbs_gen_environment_default(&env);
    env.mutation_scale_q8 = round ? (uint16_t)4096u : (uint16_t)0u;   /* x16 or none */
    bd.env = &env;
    living0 = fbs_gen_living(c->g); records0 = fbs_gen_record_count(c->g);
    CHECK_OK(fbs_gen_brood(c->g, &bd, un, written, kids, &rep));
    CHECK(rep.born + rep.culled == written && rep.culled == 0u, "born %u culled %u", rep.born, rep.culled);
    CHECK(fbs_gen_living(c->g) == living0 + rep.born, "living growth");
    CHECK(fbs_gen_record_count(c->g) == records0 + rep.born, "record growth");
    mism = 0;
    for (i = 0; i < written; i++) {
      fbs_gen_info info;
      CHECK_OK(fbs_gen_get(c->g, kids[i], &info));
      CHECK(info.sire == un[i].a && info.dam == un[i].b, "parents recorded");
      fbs_gen_genome(c->g, un[i].a, ga, 96);
      fbs_gen_genome(c->g, un[i].b, gb, 96);
      fbs_gen_genome(c->g, kids[i], gk, 96);
      for (l = 0; l < L; l++) {
        if (gk[l] != ga[l] && gk[l] != ga[L + l]) mism++;
        if (gk[L + l] != gb[l] && gk[L + l] != gb[L + l]) mism++;
      }
    }
    if (round == 0u) CHECK(mism == 0u && rep.mutations == 0u, "no mutation: %u foreign alleles, %u counted", mism, rep.mutations);
    else CHECK(mism <= rep.mutations && rep.mutations > 0u, "foreign alleles %u <= mutations %u", mism, rep.mutations);
  }
  /* retire half, then records = living + remembered dead */
  for (i = 0; i < 100u; i++) CHECK_OK(fbs_gen_retire(c->g, ids[i]));
  CHECK_ERR(fbs_gen_records(c->g, NULL, 0, &cnt), FBS_GEN_E_TRUNCATED);
  CHECK(cnt == fbs_gen_living(c->g) + 100u, "records %u = living %u + 100 dead", (unsigned)cnt, fbs_gen_living(c->g));
  free(c->mem);
  free(c);
}

/* One scripted multi-season run; returns the final snapshot (caller frees). */
static uint8_t *scripted_run(uint64_t seed, size_t *n) {
  gob_ctx *c = (gob_ctx *)malloc(sizeof *c);
  fbs_gen_id ids[300], kids[300];
  fbs_gen_contributor con[64];
  fbs_gen_union un[300];
  fbs_gen_scored top[64];
  fbs_gen_pool_desc pd;
  fbs_gen_brood_desc bd;
  fbs_gen_index ix;
  fbs_gen_environment env;
  uint32_t s, i, written;
  size_t cnt;
  uint8_t *img = NULL;
  *n = 0;
  if (!c || !gob_open(c, 1024, 300)) { free(c); return NULL; }
  fbs_gen_add_founders(c->g, 300, 0, seed, ids);
  memset(&ix, 0, sizeof ix);
  ix.weight[GOB_SIZE] = 1000; ix.weight[GOB_STRENGTH] = 1000; ix.standardize = 1;
  for (s = 1; s <= 5u; s++) {
    const fbs_gen_id *col; const int32_t *vals; uint32_t living;
    fbs_gen_set_season(c->g, s);
    fbs_gen_rank(c->g, &ix, NULL, 0, top, 40, &cnt);
    for (i = 0; i < 40u; i++) { con[i].id = top[i].id; con[i].weight = 1000u + i; }
    memset(&pd, 0, sizeof pd);
    pd.children = 150; pd.seed = seed + s; pd.kin_max_q16 = 4096; pd.redraws = 3;
    pd.assort_k = s == 3u ? 5u : 0u; pd.assort_index = &ix;
    fbs_gen_plan_pool(c->g, &pd, con, 40, un, 300, &written);
    memset(&bd, 0, sizeof bd);
    fbs_gen_environment_default(&env);
    env.mutation_scale_q8 = 1024; env.shift[GOB_STAMINA] = 80;
    env.viability_count = 1; env.viability[0].trait = GOB_DAMP;
    env.viability[0].low = 150; env.viability[0].high = 350; env.viability[0].floor_q16 = 39322;
    bd.seed = seed * 31u + s; bd.env = &env; bd.line_policy = FBS_GEN_LINE_HIGHER; bd.higher_index = &ix;
    fbs_gen_brood(c->g, &bd, un, written, kids, NULL);
    fbs_gen_column(c->g, GOB_SIZE, &col, &vals, &living);
    for (i = 0; i < living; i += 3u) if (col[i] % 5u == s % 5u) { fbs_gen_retire(c->g, col[i]); break; }
    for (i = 0; i < 30u && fbs_gen_living(c->g) > 600u; i++) {
      fbs_gen_column(c->g, GOB_SIZE, &col, &vals, &living);
      fbs_gen_retire(c->g, col[0]);
    }
  }
  img = snapshot(c->g, n);
  free(c->mem);
  free(c);
  return img;
}

static void test_determinism(void) {
  gob_ctx *c = (gob_ctx *)malloc(sizeof *c);
  fbs_gen_id ids[120], one[200], chunk[200];
  fbs_gen_contributor con[120];
  fbs_gen_union un[200];
  fbs_gen_pool_desc pd;
  fbs_gen_brood_desc bd;
  fbs_gen_brood_report rep, prep;
  fbs_gen_trait_stats pv[GOB_TRAITS], cm;
  fbs_gen_environment env;
  uint32_t i, t, written, seed;
  uint8_t *img0, *img1, *base;
  size_t n0, n1, nb;
  int8_t g1[96], g2[96];
  int32_t z1[GOB_TRAITS], z2[GOB_TRAITS];
  static const uint32_t order[4] = { 3, 1, 0, 2 };

  for (seed = 1; seed <= 3u; seed++) {
    img0 = scripted_run(seed, &n0);
    img1 = scripted_run(seed, &n1);
    CHECK(img0 && img1 && n0 == n1 && memcmp(img0, img1, n0) == 0, "seed %u: identical runs, identical bytes", seed);
    free(img1);
    img1 = scripted_run(seed + 100u, &n1);
    CHECK(img1 && (n0 != n1 || memcmp(img0, img1, n0) != 0), "different seeds differ");
    free(img0); free(img1);
  }

  if (!c || !gob_open(c, 2048, 512)) { CHECK(0, "goblin ctx"); free(c); return; }
  CHECK_OK(fbs_gen_add_founders(c->g, 120, 0, 5, ids));
  for (i = 0; i < 120u; i++) { con[i].id = ids[i]; con[i].weight = 1000u; }
  memset(&pd, 0, sizeof pd);
  pd.children = 200; pd.seed = 8; pd.kin_max_q16 = 65536; pd.redraws = 16;
  CHECK_OK(fbs_gen_plan_pool(c->g, &pd, con, 120, un, 200, &written));
  CHECK(written == 200u, "planned %u", written);
  memset(&bd, 0, sizeof bd);
  fbs_gen_environment_default(&env);
  env.mutation_scale_q8 = 2048;
  env.viability_count = 1; env.viability[0].trait = GOB_SIZE;
  env.viability[0].low = 900; env.viability[0].high = 1300; env.viability[0].floor_q16 = 16384;
  bd.seed = 42; bd.env = &env;
  base = snapshot(c->g, &nb);

  /* preview equals commit */
  CHECK_OK(fbs_gen_brood_preview(c->g, &bd, un, written, pv, GOB_TRAITS, &prep));
  CHECK(same_snapshot(c->g, base, nb), "preview changed state");
  CHECK(prep.culled > 0u && prep.born > 0u, "viability culls some: born %u culled %u", prep.born, prep.culled);
  CHECK_OK(fbs_gen_brood(c->g, &bd, un, written, one, &rep));
  CHECK(rep.born == prep.born && rep.culled == prep.culled && rep.mutations == prep.mutations &&
        rep.mean_f_q30 == prep.mean_f_q30, "preview report equals commit");
  {
    fbs_gen_id born[200];
    uint32_t k = 0;
    for (i = 0; i < written; i++) if (one[i]) born[k++] = one[i];
    CHECK(k == rep.born, "born ids");
    for (t = 0; t < GOB_TRAITS; t++) {
      CHECK_OK(fbs_gen_stats(c->g, t, born, k, &cm));
      CHECK(memcmp(&cm, &pv[t], sizeof cm) == 0, "trait %u preview stats equal commit (mean %d vs %d)", t, pv[t].mean, cm.mean);
    }
  }

  /* the same plan committed in four chunks in a scrambled order gives the same children */
  CHECK_OK(fbs_gen_load(c->g, base, nb, 0));
  for (i = 0; i < 4u; i++) {
    uint32_t k = order[i], at = k * 50u;
    CHECK_OK(fbs_gen_brood(c->g, &bd, un + at, 50, chunk + at, NULL));
  }
  for (i = 0; i < written; i++) CHECK((one[i] == 0u) == (chunk[i] == 0u), "union %u culled the same", i);
  {
    /* the one-shot children were rolled back by the load: replay the single commit on a twin
     * and compare genomes and traits per union */
    gob_ctx *d = (gob_ctx *)malloc(sizeof *d);
    if (d && gob_open(d, 2048, 512)) {
      fbs_gen_id twin[200];
      CHECK_OK(fbs_gen_load(d->g, base, nb, 0));
      CHECK_OK(fbs_gen_brood(d->g, &bd, un, written, twin, NULL));
      for (i = 0; i < written; i++) {
        if (!twin[i] || !chunk[i]) continue;
        fbs_gen_genome(d->g, twin[i], g1, 96); fbs_gen_genome(c->g, chunk[i], g2, 96);
        fbs_gen_traits(d->g, twin[i], z1, GOB_TRAITS); fbs_gen_traits(c->g, chunk[i], z2, GOB_TRAITS);
        CHECK(memcmp(g1, g2, 96) == 0 && memcmp(z1, z2, sizeof z1) == 0, "union %u child identical across chunking", i);
      }
      free(d->mem);
    }
    free(d);
  }
  free(base);
  free(c->mem);
  free(c);
}

/* ---- snapshots ----------------------------------------------------------- */

static void test_snapshots(void) {
  tiny_schema ts, tuned, renamed;
  fbs_gen_config cfg;
  fbs_gen *g, *h;
  void *mem, *mem2;
  fbs_gen_id ids[12], kids[8];
  fbs_gen_union un[8];
  fbs_gen_brood_desc bd;
  uint8_t *img, *img2, *bad;
  size_t n, n2, i;
  uint32_t k, rejected = 0, living, records;
  int32_t z0[1], z1[1];
  size_t lb, rb;

  tiny_build(&ts, 3, 120, 40);
  cfg = small_config(32, 48); cfg.record_seasons = 1;
  g = make_ctx(&cfg, &ts.schema, &mem);
  h = make_ctx(&cfg, &ts.schema, &mem2);
  CHECK(g && h, "snapshot ctxs");
  if (!g || !h) return;
  CHECK_OK(fbs_gen_add_founders(g, 12, 7, 21, ids));
  fbs_gen_set_season(g, 2);
  memset(&bd, 0, sizeof bd); bd.seed = 3; bd.line_policy = FBS_GEN_LINE_FIRST;
  for (k = 0; k < 8u; k++) { un[k].a = ids[k]; un[k].b = ids[(k + 3u) % 12u]; un[k].stream = k; un[k].kinship_q30 = 0; }
  CHECK_OK(fbs_gen_brood(g, &bd, un, 8, kids, NULL));
  CHECK_OK(fbs_gen_retire(g, ids[0]));
  CHECK_OK(fbs_gen_retire(g, ids[5]));
  CHECK_OK(fbs_gen_retire(g, kids[2]));
  CHECK_OK(fbs_gen_set_line(g, kids[3], 99));

  img = snapshot(g, &n);
  CHECK(img && n == fbs_gen_snapshot_size(g), "snapshot size");
  CHECK_ERR(fbs_gen_save(g, img, n - 1u, &n2), FBS_GEN_E_TRUNCATED);
  CHECK(n2 == n, "save reports the size");
  CHECK_OK(fbs_gen_load(h, img, n, 0));
  img2 = snapshot(h, &n2);
  CHECK(img2 && n2 == n && memcmp(img, img2, n) == 0, "save/load/save identity");
  CHECK(fbs_gen_living(h) == fbs_gen_living(g) && fbs_gen_season(h) == 2u, "loaded counts");
  CHECK_OK(fbs_gen_traits(g, kids[4], z0, 1)); CHECK_OK(fbs_gen_traits(h, kids[4], z1, 1));
  CHECK(z0[0] == z1[0], "phenotype recomputed identically");
  free(img2);

  /* every single-byte corruption is rejected and leaves the target unchanged */
  bad = (uint8_t *)malloc(n);
  img2 = snapshot(h, &n2);
  for (i = 0; i < n; i++) {
    memcpy(bad, img, n);
    bad[i] ^= (uint8_t)(1u + (i % 255u));
    if (fbs_gen_load(h, bad, n, FBS_GEN_LOAD_RETUNE) != FBS_GEN_OK) rejected++;
  }
  CHECK(rejected == n, "%u of %u single-byte corruptions rejected", rejected, (unsigned)n);
  CHECK(same_snapshot(h, img2, n2), "rejected loads left state unchanged");

  /* semantic corruptions with a valid checksum */
  living = get32(img + 52); records = get32(img + 56);
  lb = 4u + 4u * 1u + 2u * 3u; rb = 36u + 4u;
  {
    uint8_t *recs = img + 64u + living * lb, *liv = img + 64u;
    /* a child older than its sire */
    memcpy(bad, img, n); put32(bad + (recs - img) + (records - 1u) * rb + 4u, get32(recs + (records - 1u) * rb)); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* duplicate serials in records */
    memcpy(bad, img, n); put32(bad + (recs - img) + rb, get32(recs)); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* a living entry whose record is dead */
    memcpy(bad, img, n); put32(bad + (liv - img), ids[0]); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* duplicate living serials */
    memcpy(bad, img, n); put32(bad + (liv - img) + lb, get32(liv)); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* allele out of range */
    memcpy(bad, img, n); bad[(liv - img) + 8u] = 9; reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* a record born after the season */
    memcpy(bad, img, n); put32(bad + (recs - img) + 12u, 3u); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* next serial not above every serial */
    memcpy(bad, img, n); put32(bad + 48u, get32(recs + (records - 1u) * rb)); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* unknown flag bits */
    memcpy(bad, img, n); put32(bad + (recs - img) + 24u, 0x80u); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* F above 1 */
    memcpy(bad, img, n); put32(bad + (recs - img) + 20u, Q30 + 1u); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_CORRUPT);
    /* next serial at the top: one more birth uses the last serial, then E_FULL, never a wrap */
    {
      void *mem3;
      fbs_gen *w = make_ctx(&cfg, &ts.schema, &mem3);
      fbs_gen_id last = 0, more = 0;
      memcpy(bad, img, n); put32(bad + 48u, 0xFFFFFFFFu); reseal(bad, n);
      if (w) {
        uint8_t *wimg;
        size_t wn;
        CHECK_OK(fbs_gen_load(w, bad, n, 0));
        CHECK_ERR(fbs_gen_add_founders(w, 2, 0, 1, NULL), FBS_GEN_E_FULL);
        CHECK_OK(fbs_gen_add_founders(w, 1, 0, 1, &last));
        CHECK(last == 0xFFFFFFFFu, "last serial %u", last);
        wimg = snapshot(w, &wn);
        CHECK_ERR(fbs_gen_add_founders(w, 1, 0, 1, &more), FBS_GEN_E_FULL);
        CHECK(wimg && same_snapshot(w, wimg, wn), "exhausted serials refuse without change");
        CHECK(wimg && fbs_gen_load(w, wimg, wn, 0) == FBS_GEN_OK, "an exhausted context reloads");
        free(wimg);
        free(mem3);
      }
    }
    /* snapshot version */
    memcpy(bad, img, n); put32(bad + 8u, 2u); reseal(bad, n);
    CHECK_ERR(fbs_gen_load(h, bad, n, 0), FBS_GEN_E_SCHEMA);
  }
  CHECK(same_snapshot(h, img2, n2), "semantic rejections left state unchanged");
  free(img2);
  free(bad);
  free(mem2);

  /* retune: a different weight is accepted only with FBS_GEN_LOAD_RETUNE */
  tuned = ts;
  tuned.schema.traits = &tuned.trait; tuned.schema.loci = tuned.loci; tuned.schema.effects = tuned.effects;
  for (k = 0; k < 3u; k++) { tuned.loci[k].key = tuned.keys[k]; tuned.effects[k].weight = 240; }
  h = make_ctx(&cfg, &tuned.schema, &mem2);
  CHECK(h != NULL, "tuned ctx");
  if (h) {
    uint64_t s1, t1, s2, t2;
    fbs_gen_schema_hashes(g, &s1, &t1); fbs_gen_schema_hashes(h, &s2, &t2);
    CHECK(s1 == s2 && t1 != t2, "tuning changes only the tuning hash");
    CHECK_ERR(fbs_gen_load(h, img, n, 0), FBS_GEN_E_SCHEMA);
    CHECK(fbs_gen_living(h) == 0u, "refused retune left the context empty");
    CHECK_OK(fbs_gen_load(h, img, n, FBS_GEN_LOAD_RETUNE));
    {
      int32_t a0[1], a1[1];
      CHECK_OK(fbs_gen_additive(g, kids[4], a0, 1)); CHECK_OK(fbs_gen_additive(h, kids[4], a1, 1));
      CHECK(a1[0] - 10000 == 2 * (a0[0] - 10000), "retuned additive value doubles: %d vs %d", a0[0], a1[0]);
    }
    free(mem2);
  }
  renamed = ts;
  renamed.schema.traits = &renamed.trait; renamed.schema.loci = renamed.loci; renamed.schema.effects = renamed.effects;
  for (k = 0; k < 3u; k++) renamed.loci[k].key = renamed.keys[k];
  renamed.keys[2][0] = 'M';
  h = make_ctx(&cfg, &renamed.schema, &mem2);
  CHECK(h != NULL, "renamed ctx");
  if (h) {
    CHECK_ERR(fbs_gen_load(h, img, n, FBS_GEN_LOAD_RETUNE), FBS_GEN_E_SCHEMA);
    free(mem2);
  }
  /* a context with record_traits off cannot take the image */
  cfg.record_traits = 0;
  h = make_ctx(&cfg, &ts.schema, &mem2);
  if (h) {
    int32_t bz[1];
    CHECK_ERR(fbs_gen_load(h, img, n, 0), FBS_GEN_E_SCHEMA);
    CHECK_ERR(fbs_gen_birth_traits(h, ids[1], bz, 1), FBS_GEN_E_STATE);
    free(mem2);
  }
  /* smaller capacity: E_FULL */
  cfg = small_config(4, 48); cfg.record_seasons = 1;
  h = make_ctx(&cfg, &ts.schema, &mem2);
  if (h) { CHECK_ERR(fbs_gen_load(h, img, n, 0), FBS_GEN_E_FULL); free(mem2); }
  free(img);
  free(mem);
}

int main(void) {
  test_bounds();
  test_kinship();
  test_eviction();
  test_conservation();
  test_determinism();
  test_snapshots();
  printf("%d checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
