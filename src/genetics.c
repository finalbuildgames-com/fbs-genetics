/*
 * src/genetics.c — FinalBuildSystems heritable traits. Implements include/fbs/genetics.h.
 *
 * Original work. The model follows textbook infinitesimal/continuum-of-alleles quantitative
 * genetics (Falconer & Mackay; Lynch & Walsh; Wright's path kinship; Haldane's map function for
 * the authoring tool). The PRNG is public domain: Blackman & Vigna's xoshiro128** 1.1 seeded by
 * Vigna's splitmix64.
 *
 * C99, <string.h> only. No allocation, no libm, no globals, no threads. Doubles appear only in
 * outputs (variances for forecasts and standardized scores) and never feed back into state; the
 * module is built with -ffp-contract=off so the compiler does not fuse those operations. The
 * tests check that identical call sequences give identical snapshot bytes on one build; no
 * shipped test compares native and WASM results.
 *
 * Internal contracts the header only implies:
 *
 * STREAMS. The child of union u in a brood with seed s draws, in this order: the first gamete's
 * starting strand, one 32-bit draw per later locus for crossover, a mutation draw per locus whose
 * effective rate is non-zero (plus the step and sign draws when it fires), then the second gamete
 * likewise, then two draws per trait with non-zero noise (Irwin-Hall of four 16-bit uniforms), then
 * one draw per viability rule. The stream is xoshiro128**(splitmix64-seeded from
 * s ^ splitmix64(u.stream)).
 *
 * KINSHIP. f(a, a) = (1 + F_a) / 2; for serial b > a, f(a, b) = (f(a, sire_b) + f(a, dam_b)) / 2.
 * Each side carries its generation distance from the query; a node at distance pedigree_depth is
 * not expanded (its parents count as unknown). Sums are exact in uint64 and halved with a
 * truncating shift. A query lays both ancestor trees out in heap order (position p has parents
 * 2p+1 and 2p+2, generation floor(log2(p+1))) and memoizes the recursion densely over (position in
 * a's tree, position in b's tree), since expanding the younger node always stays inside its own
 * tree. f(x, y) can only be non-zero when the two subtrees share a serial, so each node carries a
 * bit mask of the shared serials below it (up to 64; more disables the pruning) and disjoint
 * pairs return 0 without recursion. The persistent memo caches top-level results by (a, b) only; it is cleared whenever a
 * record disappears (eviction, forget, load), so neither cache ever changes a result.
 *
 * RECORDS. Records are kept sorted by serial (serials are birth order and records are appended),
 * looked up by binary search, and compacted in place on eviction and forget. Each record caches
 * its parents' record indices so pedigree walks need no search; compaction remaps them.
 */
#include "fbs/genetics.h"

#include <string.h>

/* ------------------------------------------------------------------------- */
/* Internal types                                                            */
/* ------------------------------------------------------------------------- */

#define GEN_MAGIC 0x46424731u /* "FBG1" */
#define GEN_NONE  0xFFFFFFFFu
#define GEN_HEADER_BYTES 64u
#define GEN_RECORD_BYTES 36u
#define GEN_ANC_MAX 511u /* heap positions for depth 8: 2^9 - 1 */

typedef struct gen_trait {
  int32_t base, lo, hi, env_sd, depression;
  uint32_t key_off, key_len;
} gen_trait;

typedef struct gen_locus {
  int8_t amin, amax;
  uint8_t step, fcount;
  uint16_t h, rec, mut, pad;
  int8_t fa[8];
  uint16_t fw[8];
  uint32_t fw_total;
  uint32_t key_off, key_len;
  uint32_t eff_begin, eff_end;
} gen_locus;

typedef struct gen_effect { uint16_t locus; uint8_t trait; int8_t match; int32_t weight; } gen_effect;

typedef struct gen_record {
  uint32_t serial, sire, dam, born, line, f, flags, children, honours, slot;
  uint32_t sire_r, dam_r;        /* record indices of the parents, GEN_NONE when not recorded */
} gen_record;

typedef struct gen_memo { uint32_t a, b, g, v; } gen_memo;

typedef struct gen_rng { uint32_t s[4]; } gen_rng;

struct fbs_gen {
  uint32_t magic;
  fbs_gen_config cfg;
  uint32_t L, T, E, stride;
  uint64_t hash_structure, hash_tuning;
  gen_trait *traits;
  gen_locus *loci;
  gen_effect *effects;
  char *keys;
  uint32_t season, next_serial;
  uint32_t living, records;
  uint32_t brood_counter;
  /* living, dense slots */
  uint32_t *l_serial, *l_f, *l_stamp;
  int8_t *l_al;                  /* stride bytes per slot */
  int32_t *l_z, *l_a, *l_e;      /* SoA: [t * max_living + slot] */
  /* records */
  gen_record *rec;
  int32_t *rec_z;                /* [r * T + t] when record_traits */
  uint32_t *rec_map;             /* compaction scratch: old record index -> new */
  /* kinship memo */
  gen_memo *memo;
  /* brood scratch */
  int8_t *b_al;
  int32_t *b_e, *b_z, *b_a;      /* [i * T + t] */
  uint32_t *b_f, *b_line;
  uint8_t *b_ok;
  /* plan scratch */
  uint64_t *p_cum;
  int64_t *p_score;
  uint32_t *p_order, *p_pos, *p_nb, *p_slot;
  /* rank / load scratch */
  fbs_gen_scored *r_buf;
  int64_t *r_score;
  uint32_t *r_slot;
  /* kinship scratch: two heap-ordered ancestor trees and a stamped dense memo [H * H] */
  uint32_t k_heap;               /* H = 2^(depth+1) - 1 */
  uint32_t k_stamp;
  uint64_t *k_memo;
  uint32_t k_ser[2][GEN_ANC_MAX], k_f[2][GEN_ANC_MAX], k_rec[2][GEN_ANC_MAX];
  uint8_t k_has[2][GEN_ANC_MAX];
  uint64_t k_mask[2][GEN_ANC_MAX];   /* shared ancestors (by index into k_shared) in each subtree */
  uint32_t k_shared[64];
  uint32_t anc_a[GEN_ANC_MAX], anc_b[GEN_ANC_MAX];
};

typedef struct gen_layout {
  size_t traits, loci, effects, keys;
  size_t l_serial, l_f, l_stamp, l_al, l_z, l_a, l_e;
  size_t rec, rec_z, rec_map, memo, k_memo;
  size_t b_al, b_e, b_z, b_a, b_f, b_line, b_ok;
  size_t p_cum, p_score, p_order, p_pos, p_nb, p_slot;
  size_t r_buf, r_score, r_slot;
  size_t total;
} gen_layout;

/* ------------------------------------------------------------------------- */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------- */

const char *fbs_gen_status_name(int status) {
  switch (status) {
    case FBS_GEN_OK: return "FBS_GEN_OK";
    case FBS_GEN_E_INVALID: return "FBS_GEN_E_INVALID";
    case FBS_GEN_E_RANGE: return "FBS_GEN_E_RANGE";
    case FBS_GEN_E_FULL: return "FBS_GEN_E_FULL";
    case FBS_GEN_E_NOT_FOUND: return "FBS_GEN_E_NOT_FOUND";
    case FBS_GEN_E_DEAD: return "FBS_GEN_E_DEAD";
    case FBS_GEN_E_STATE: return "FBS_GEN_E_STATE";
    case FBS_GEN_E_SCHEMA: return "FBS_GEN_E_SCHEMA";
    case FBS_GEN_E_TRUNCATED: return "FBS_GEN_E_TRUNCATED";
    case FBS_GEN_E_MEMORY: return "FBS_GEN_E_MEMORY";
    case FBS_GEN_E_CORRUPT: return "FBS_GEN_E_CORRUPT";
    default: return "FBS_GEN_E_UNKNOWN";
  }
}

unsigned fbs_gen_version(void) { return FBS_GEN_VERSION; }

static int ctx_ok(const fbs_gen *c) { return c != NULL && c->magic == GEN_MAGIC; }

static int64_t round_div(int64_t num, int64_t den) {
  /* den > 0; half away from zero */
  if (num >= 0) return (num + den / 2) / den;
  return -((-num + den / 2) / den);
}

static int64_t round_d(double v) {
  if (v >= 9.2e18) return INT64_MAX;
  if (v <= -9.2e18) return INT64_MIN;
  return v >= 0.0 ? (int64_t)(v + 0.5) : -(int64_t)(-v + 0.5);
}

static int32_t clamp32(int64_t v, int32_t lo, int32_t hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return (int32_t)v;
}

static int32_t sat32(int64_t v) { return clamp32(v, INT32_MIN, INT32_MAX); }

static uint64_t isqrt_u64(uint64_t v) {
  uint64_t r = 0, bit = (uint64_t)1 << 62;
  while (bit > v) bit >>= 2;
  while (bit) {
    if (v >= r + bit) { v -= r + bit; r = (r >> 1) + bit; }
    else r >>= 1;
    bit >>= 2;
  }
  return r;
}

static int32_t sd_from_var(double var) {
  uint64_t v;
  if (!(var > 0.0)) return 0;
  if (var > 4.0e18) var = 4.0e18;
  v = (uint64_t)(var + 0.5);
  return sat32((int64_t)isqrt_u64(v));
}

static uint64_t fnv_bytes(uint64_t h, const void *p, size_t n) {
  const uint8_t *b = (const uint8_t *)p;
  size_t i;
  for (i = 0; i < n; i++) { h ^= b[i]; h *= 0x100000001B3ULL; }
  return h;
}
static uint64_t fnv_u32(uint64_t h, uint32_t v) {
  uint8_t b[4];
  b[0] = (uint8_t)v; b[1] = (uint8_t)(v >> 8); b[2] = (uint8_t)(v >> 16); b[3] = (uint8_t)(v >> 24);
  return fnv_bytes(h, b, 4);
}
#define FNV_INIT 0xCBF29CE484222325ULL

static void put_u32(uint8_t *p, uint32_t v) {
  p[0] = (uint8_t)v; p[1] = (uint8_t)(v >> 8); p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
}
static void put_u64(uint8_t *p, uint64_t v) { put_u32(p, (uint32_t)v); put_u32(p + 4, (uint32_t)(v >> 32)); }
static uint32_t get_u32(const uint8_t *p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static uint64_t get_u64(const uint8_t *p) { return (uint64_t)get_u32(p) | ((uint64_t)get_u32(p + 4) << 32); }

/* ------------------------------------------------------------------------- */
/* PRNG: xoshiro128** 1.1 seeded by splitmix64                              */
/* ------------------------------------------------------------------------- */

static uint32_t rotl32(uint32_t x, unsigned k) { return (uint32_t)((x << k) | (x >> (32u - k))); }

static uint64_t splitmix64(uint64_t *state) {
  uint64_t z;
  *state += 0x9E3779B97F4A7C15ULL;
  z = *state;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static void rng_seed(gen_rng *r, uint64_t seed) {
  uint64_t state = seed, w;
  w = splitmix64(&state);
  r->s[0] = (uint32_t)(w & 0xFFFFFFFFu); r->s[1] = (uint32_t)(w >> 32);
  w = splitmix64(&state);
  r->s[2] = (uint32_t)(w & 0xFFFFFFFFu); r->s[3] = (uint32_t)(w >> 32);
}

static void rng_stream(gen_rng *r, uint64_t seed, uint32_t stream) {
  uint64_t st = stream;
  rng_seed(r, seed ^ splitmix64(&st));
}

static uint32_t rng_next(gen_rng *r) {
  uint32_t result = rotl32((uint32_t)(r->s[1] * 5u), 7u) * 9u;
  uint32_t t = (uint32_t)(r->s[1] << 9);
  r->s[2] ^= r->s[0];
  r->s[3] ^= r->s[1];
  r->s[1] ^= r->s[2];
  r->s[0] ^= r->s[3];
  r->s[2] ^= t;
  r->s[3] = rotl32(r->s[3], 11u);
  return result;
}

static uint32_t rng_below(gen_rng *r, uint32_t bound) {
  uint32_t threshold, x;
  if (bound == 0u) return 0u;
  threshold = (uint32_t)(0u - bound) % bound;
  do { x = rng_next(r); } while (x < threshold);
  return x % bound;
}

static uint64_t rng_below64(gen_rng *r, uint64_t bound) {
  uint64_t threshold, x;
  if (bound == 0u) return 0u;
  if (bound <= 0xFFFFFFFFull) return rng_below(r, (uint32_t)bound);
  threshold = (0u - bound) % bound;
  do {
    x = (uint64_t)rng_next(r) << 32;
    x |= rng_next(r);
  } while (x < threshold);
  return x % bound;
}

/* Irwin-Hall(4) scaled to SD sd (milli-units); bounded at +-2*sqrt(3)*sd. */
static int32_t rng_noise(gen_rng *r, int64_t sd) {
  uint32_t x, y;
  int64_t s;
  if (sd <= 0) return 0;
  x = rng_next(r); y = rng_next(r);
  s = (int64_t)(x & 0xFFFFu) + (int64_t)(x >> 16) + (int64_t)(y & 0xFFFFu) + (int64_t)(y >> 16) - 131070;
  return (int32_t)round_div(s * 113512 * sd, (int64_t)1 << 32); /* 113512 = sqrt(3) in Q16 */
}

/* ------------------------------------------------------------------------- */
/* Schema validation and layout                                              */
/* ------------------------------------------------------------------------- */

static int32_t effect_value(const gen_effect *e, int32_t a) {
  if (e->match == FBS_GEN_ANY) return e->weight * a; /* validated to fit */
  return a == e->match ? e->weight : 0;
}

fbs_gen_status fbs_gen_schema_check(const fbs_gen_schema *s) {
  uint32_t i, j;
  int64_t worst[FBS_GEN_MAX_TRAITS];
  if (!s || !s->traits || !s->loci) return FBS_GEN_E_INVALID;
  if (s->trait_count < 1u || s->trait_count > FBS_GEN_MAX_TRAITS) return FBS_GEN_E_SCHEMA;
  if (s->locus_count < 1u || s->locus_count > FBS_GEN_MAX_LOCI) return FBS_GEN_E_SCHEMA;
  if (s->effect_count > FBS_GEN_MAX_EFFECTS) return FBS_GEN_E_SCHEMA;
  if (s->effect_count > 0u && !s->effects) return FBS_GEN_E_INVALID;
  for (i = 0; i < s->trait_count; i++) {
    const fbs_gen_trait_desc *t = &s->traits[i];
    if (!t->key || t->key_len < 1u || t->key_len > FBS_GEN_KEY_MAX) return FBS_GEN_E_SCHEMA;
    for (j = 0; j < i; j++)
      if (s->traits[j].key_len == t->key_len && memcmp(s->traits[j].key, t->key, t->key_len) == 0)
        return FBS_GEN_E_SCHEMA;
    if (t->lo > t->base || t->base > t->hi) return FBS_GEN_E_SCHEMA;
    if (t->env_sd < 0 || t->env_sd > FBS_GEN_ENV_SD_MAX) return FBS_GEN_E_RANGE;
    if (t->depression > (1 << 30) || t->depression < -(1 << 30)) return FBS_GEN_E_RANGE;
    worst[i] = t->base < 0 ? -(int64_t)t->base : (int64_t)t->base;
  }
  for (i = 0; i < s->locus_count; i++) {
    const fbs_gen_locus_desc *l = &s->loci[i];
    if (!l->key || l->key_len < 1u || l->key_len > FBS_GEN_KEY_MAX) return FBS_GEN_E_SCHEMA;
    for (j = 0; j < i; j++)
      if (s->loci[j].key_len == l->key_len && memcmp(s->loci[j].key, l->key, l->key_len) == 0)
        return FBS_GEN_E_SCHEMA;
    if (l->amin <= FBS_GEN_ANY || l->amin > l->amax) return FBS_GEN_E_SCHEMA;
    if (l->h > 256u) return FBS_GEN_E_SCHEMA;
    if (l->recombination_q16 > 32768u) return FBS_GEN_E_SCHEMA;
    if (i == 0 && l->recombination_q16 != 32768u) return FBS_GEN_E_SCHEMA;
    if (l->mutation_q16 > 0u && (l->mutation_step < 1u || l->mutation_step > 127u)) return FBS_GEN_E_SCHEMA;
    if (l->founder_count < 1u || l->founder_count > 8u) return FBS_GEN_E_SCHEMA;
    for (j = 0; j < l->founder_count; j++) {
      if (l->founder[j].weight == 0u) return FBS_GEN_E_SCHEMA;
      if (l->founder[j].allele < l->amin || l->founder[j].allele > l->amax) return FBS_GEN_E_SCHEMA;
    }
  }
  for (i = 0; i < s->effect_count; i++) {
    const fbs_gen_effect *e = &s->effects[i];
    const fbs_gen_locus_desc *l;
    int64_t w, m;
    if (e->locus >= s->locus_count || e->trait >= s->trait_count) return FBS_GEN_E_SCHEMA;
    l = &s->loci[e->locus];
    if (e->match != FBS_GEN_ANY && (e->match < l->amin || e->match > l->amax)) return FBS_GEN_E_SCHEMA;
    w = e->weight < 0 ? -(int64_t)e->weight : (int64_t)e->weight;
    if (e->match == FBS_GEN_ANY) {
      int64_t a = -(int64_t)l->amin > (int64_t)l->amax ? -(int64_t)l->amin : (int64_t)l->amax;
      if (a < 0) a = -a;
      m = w * a;
    } else {
      m = w;
    }
    if (m > (1 << 30)) return FBS_GEN_E_RANGE;
    worst[e->trait] += 2 * m;
    if (worst[e->trait] > (1 << 30)) return FBS_GEN_E_RANGE;
  }
  return FBS_GEN_OK;
}

static int config_resolve(const fbs_gen_config *in, fbs_gen_config *out) {
  if (!in) return 0;
  *out = *in;
  if (out->pedigree_depth == 0u) out->pedigree_depth = 5u;
  if (out->record_seasons == 0u) out->record_seasons = 8u;
  if (out->memo_entries == 0u) out->memo_entries = 4096u;
  if (out->max_living < 1u || out->max_living > FBS_GEN_MAX_LIVING) return 0;
  if (out->max_records < out->max_living || out->max_records > FBS_GEN_MAX_RECORDS) return 0;
  if (out->pedigree_depth > FBS_GEN_MAX_DEPTH) return 0;
  if (out->memo_entries > (1u << 24) || (out->memo_entries & (out->memo_entries - 1u)) != 0u) return 0;
  if (out->max_brood < 1u || out->max_brood > FBS_GEN_MAX_BROOD) return 0;
  if (out->max_contributors < 1u || out->max_contributors > out->max_living) return 0;
  if (out->record_traits > 1u) return 0;
  return 1;
}

/* Offsets are computed in uint64_t so a 32-bit size_t cannot wrap; layout_build returns 0 when
 * the total does not fit size_t, and then no offset is used. */
static size_t layout_take(uint64_t *at, uint64_t bytes) {
  uint64_t o = (*at + 7u) & ~(uint64_t)7u;
  *at = o + bytes;
  return (size_t)o;
}

static size_t layout_build(const fbs_gen_config *cfg, const fbs_gen_schema *s, gen_layout *lay) {
  uint64_t at = ((uint64_t)sizeof(fbs_gen) + 7u) & ~(uint64_t)7u;
  uint64_t keys = 0, L = s->locus_count, T = s->trait_count, stride, cap = cfg->max_living;
  uint64_t rc = cfg->max_records, bc = cfg->max_brood, pc = cfg->max_contributors;
  uint32_t i;
  for (i = 0; i < s->trait_count; i++) keys += s->traits[i].key_len;
  for (i = 0; i < s->locus_count; i++) keys += s->loci[i].key_len;
  stride = (2u * L + 15u) & ~(uint64_t)15u;
  lay->traits = layout_take(&at, sizeof(gen_trait) * T);
  lay->loci = layout_take(&at, sizeof(gen_locus) * L);
  lay->effects = layout_take(&at, sizeof(gen_effect) * (s->effect_count ? s->effect_count : 1u));
  lay->keys = layout_take(&at, keys);
  lay->l_serial = layout_take(&at, 4u * cap);
  lay->l_f = layout_take(&at, 4u * cap);
  lay->l_stamp = layout_take(&at, 4u * cap);
  lay->l_al = layout_take(&at, stride * cap);
  lay->l_z = layout_take(&at, 4u * T * cap);
  lay->l_a = layout_take(&at, 4u * T * cap);
  lay->l_e = layout_take(&at, 4u * T * cap);
  lay->rec = layout_take(&at, sizeof(gen_record) * rc);
  lay->rec_z = layout_take(&at, cfg->record_traits ? 4u * T * rc : 0u);
  lay->rec_map = layout_take(&at, 4u * rc);
  lay->memo = layout_take(&at, sizeof(gen_memo) * cfg->memo_entries);
  {
    uint64_t h = ((uint64_t)1 << (cfg->pedigree_depth + 1u)) - 1u;
    lay->k_memo = layout_take(&at, 8u * h * h);
  }
  lay->b_al = layout_take(&at, stride * bc);
  lay->b_e = layout_take(&at, 4u * T * bc);
  lay->b_z = layout_take(&at, 4u * T * bc);
  lay->b_a = layout_take(&at, 4u * T * bc);
  lay->b_f = layout_take(&at, 4u * bc);
  lay->b_line = layout_take(&at, 4u * bc);
  lay->b_ok = layout_take(&at, bc);
  lay->p_cum = layout_take(&at, 8u * pc);
  lay->p_score = layout_take(&at, 8u * pc);
  lay->p_order = layout_take(&at, 4u * pc);
  lay->p_pos = layout_take(&at, 4u * pc);
  lay->p_nb = layout_take(&at, 4u * pc);
  lay->p_slot = layout_take(&at, 4u * pc);
  lay->r_buf = layout_take(&at, sizeof(fbs_gen_scored) * cap);
  lay->r_score = layout_take(&at, 8u * cap);
  lay->r_slot = layout_take(&at, 4u * cap);
  at = (at + 7u) & ~(uint64_t)7u;
  if (at > (uint64_t)SIZE_MAX) { lay->total = 0; return 0; }
  lay->total = (size_t)at;
  return lay->total;
}

fbs_gen_config fbs_gen_config_default(void) {
  fbs_gen_config c;
  c.max_living = 1024u;
  c.max_records = 2048u;
  c.pedigree_depth = 5u;
  c.record_seasons = 8u;
  c.memo_entries = 4096u;
  c.max_brood = 2048u;
  c.max_contributors = 1024u;
  c.record_traits = 1u;
  return c;
}

size_t fbs_gen_memory_for(const fbs_gen_config *cfg, const fbs_gen_schema *schema) {
  fbs_gen_config r;
  gen_layout lay;
  if (!config_resolve(cfg, &r)) return 0;
  if (fbs_gen_schema_check(schema) != FBS_GEN_OK) return 0;
  return layout_build(&r, schema, &lay);
}

static void compute_hashes(fbs_gen *c) {
  uint64_t hs = FNV_INIT, ht = FNV_INIT;
  uint32_t i, j;
  hs = fnv_u32(hs, c->L); hs = fnv_u32(hs, c->T);
  for (i = 0; i < c->T; i++) {
    const gen_trait *t = &c->traits[i];
    hs = fnv_u32(hs, t->key_len); hs = fnv_bytes(hs, c->keys + t->key_off, t->key_len);
    ht = fnv_u32(ht, (uint32_t)t->base); ht = fnv_u32(ht, (uint32_t)t->lo); ht = fnv_u32(ht, (uint32_t)t->hi);
    ht = fnv_u32(ht, (uint32_t)t->env_sd); ht = fnv_u32(ht, (uint32_t)t->depression);
  }
  for (i = 0; i < c->L; i++) {
    const gen_locus *l = &c->loci[i];
    hs = fnv_u32(hs, l->key_len); hs = fnv_bytes(hs, c->keys + l->key_off, l->key_len);
    hs = fnv_u32(hs, (uint32_t)(int32_t)l->amin); hs = fnv_u32(hs, (uint32_t)(int32_t)l->amax);
    ht = fnv_u32(ht, l->h); ht = fnv_u32(ht, l->rec); ht = fnv_u32(ht, l->mut); ht = fnv_u32(ht, l->step);
    ht = fnv_u32(ht, l->fcount);
    for (j = 0; j < l->fcount; j++) { ht = fnv_u32(ht, (uint32_t)(int32_t)l->fa[j]); ht = fnv_u32(ht, l->fw[j]); }
  }
  ht = fnv_u32(ht, c->E);
  for (i = 0; i < c->E; i++) {
    const gen_effect *e = &c->effects[i];
    ht = fnv_u32(ht, e->locus); ht = fnv_u32(ht, e->trait);
    ht = fnv_u32(ht, (uint32_t)(int32_t)e->match); ht = fnv_u32(ht, (uint32_t)e->weight);
  }
  c->hash_structure = hs;
  c->hash_tuning = ht;
}

fbs_gen_status fbs_gen_create(const fbs_gen_config *cfg, const fbs_gen_schema *schema,
                              void *memory, size_t bytes, fbs_gen **out) {
  fbs_gen_config r;
  gen_layout lay;
  fbs_gen *c;
  uint8_t *base;
  fbs_gen_status st;
  uint32_t i, j, koff = 0;
  uint32_t count[FBS_GEN_MAX_LOCI + 1];
  if (!cfg || !schema || !memory || !out) return FBS_GEN_E_INVALID;
  if (((uintptr_t)memory & 7u) != 0u) return FBS_GEN_E_INVALID;
  if (!config_resolve(cfg, &r)) return FBS_GEN_E_RANGE;
  st = fbs_gen_schema_check(schema);
  if (st != FBS_GEN_OK) return st;
  if (layout_build(&r, schema, &lay) == 0u) return FBS_GEN_E_RANGE;   /* does not fit size_t */
  if (bytes < lay.total) return FBS_GEN_E_MEMORY;
  base = (uint8_t *)memory;
  memset(base, 0, lay.total);
  c = (fbs_gen *)memory;
  c->cfg = r;
  c->L = schema->locus_count; c->T = schema->trait_count; c->E = schema->effect_count;
  c->stride = (2u * c->L + 15u) & ~15u;
  c->traits = (gen_trait *)(void *)(base + lay.traits);
  c->loci = (gen_locus *)(void *)(base + lay.loci);
  c->effects = (gen_effect *)(void *)(base + lay.effects);
  c->keys = (char *)(base + lay.keys);
  c->l_serial = (uint32_t *)(void *)(base + lay.l_serial);
  c->l_f = (uint32_t *)(void *)(base + lay.l_f);
  c->l_stamp = (uint32_t *)(void *)(base + lay.l_stamp);
  c->l_al = (int8_t *)(base + lay.l_al);
  c->l_z = (int32_t *)(void *)(base + lay.l_z);
  c->l_a = (int32_t *)(void *)(base + lay.l_a);
  c->l_e = (int32_t *)(void *)(base + lay.l_e);
  c->rec = (gen_record *)(void *)(base + lay.rec);
  c->rec_z = (int32_t *)(void *)(base + lay.rec_z);
  c->rec_map = (uint32_t *)(void *)(base + lay.rec_map);
  c->memo = (gen_memo *)(void *)(base + lay.memo);
  c->k_memo = (uint64_t *)(void *)(base + lay.k_memo);
  c->k_heap = (1u << (r.pedigree_depth + 1u)) - 1u;
  c->k_stamp = 0;
  c->b_al = (int8_t *)(base + lay.b_al);
  c->b_e = (int32_t *)(void *)(base + lay.b_e);
  c->b_z = (int32_t *)(void *)(base + lay.b_z);
  c->b_a = (int32_t *)(void *)(base + lay.b_a);
  c->b_f = (uint32_t *)(void *)(base + lay.b_f);
  c->b_line = (uint32_t *)(void *)(base + lay.b_line);
  c->b_ok = base + lay.b_ok;
  c->p_cum = (uint64_t *)(void *)(base + lay.p_cum);
  c->p_score = (int64_t *)(void *)(base + lay.p_score);
  c->p_order = (uint32_t *)(void *)(base + lay.p_order);
  c->p_pos = (uint32_t *)(void *)(base + lay.p_pos);
  c->p_nb = (uint32_t *)(void *)(base + lay.p_nb);
  c->p_slot = (uint32_t *)(void *)(base + lay.p_slot);
  c->r_buf = (fbs_gen_scored *)(void *)(base + lay.r_buf);
  c->r_score = (int64_t *)(void *)(base + lay.r_score);
  c->r_slot = (uint32_t *)(void *)(base + lay.r_slot);
  for (i = 0; i < c->T; i++) {
    const fbs_gen_trait_desc *d = &schema->traits[i];
    gen_trait *t = &c->traits[i];
    t->base = d->base; t->lo = d->lo; t->hi = d->hi; t->env_sd = d->env_sd; t->depression = d->depression;
    t->key_off = koff; t->key_len = (uint32_t)d->key_len;
    memcpy(c->keys + koff, d->key, d->key_len); koff += (uint32_t)d->key_len;
  }
  for (i = 0; i < c->L; i++) {
    const fbs_gen_locus_desc *d = &schema->loci[i];
    gen_locus *l = &c->loci[i];
    l->amin = d->amin; l->amax = d->amax; l->h = d->h; l->rec = d->recombination_q16;
    l->mut = d->mutation_q16; l->step = d->mutation_q16 ? d->mutation_step : 0u;
    l->fcount = d->founder_count; l->fw_total = 0;
    for (j = 0; j < d->founder_count; j++) {
      l->fa[j] = d->founder[j].allele; l->fw[j] = d->founder[j].weight; l->fw_total += d->founder[j].weight;
    }
    l->key_off = koff; l->key_len = (uint32_t)d->key_len;
    memcpy(c->keys + koff, d->key, d->key_len); koff += (uint32_t)d->key_len;
  }
  /* Stable counting sort of effects by locus. */
  memset(count, 0, sizeof count);
  for (i = 0; i < c->E; i++) count[schema->effects[i].locus + 1u]++;
  for (i = 0; i < c->L; i++) count[i + 1u] += count[i];
  for (i = 0; i < c->L; i++) { c->loci[i].eff_begin = count[i]; c->loci[i].eff_end = count[i + 1u]; }
  for (i = 0; i < c->E; i++) {
    const fbs_gen_effect *e = &schema->effects[i];
    gen_effect *d = &c->effects[count[e->locus]++];
    d->locus = e->locus; d->trait = e->trait; d->match = e->match; d->weight = e->weight;
  }
  compute_hashes(c);
  c->season = 0; c->next_serial = 1; c->living = 0; c->records = 0; c->brood_counter = 0;
  c->magic = GEN_MAGIC;
  *out = c;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_schema_hashes(const fbs_gen *c, uint64_t *structure, uint64_t *tuning) {
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  if (structure) *structure = c->hash_structure;
  if (tuning) *tuning = c->hash_tuning;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_trait_find(const fbs_gen *c, const char *key, size_t len, uint32_t *trait) {
  uint32_t i;
  if (!ctx_ok(c) || !key || !trait) return FBS_GEN_E_INVALID;
  for (i = 0; i < c->T; i++)
    if (c->traits[i].key_len == len && memcmp(c->keys + c->traits[i].key_off, key, len) == 0) {
      *trait = i; return FBS_GEN_OK;
    }
  return FBS_GEN_E_NOT_FOUND;
}

fbs_gen_status fbs_gen_locus_find(const fbs_gen *c, const char *key, size_t len, uint32_t *locus) {
  uint32_t i;
  if (!ctx_ok(c) || !key || !locus) return FBS_GEN_E_INVALID;
  for (i = 0; i < c->L; i++)
    if (c->loci[i].key_len == len && memcmp(c->keys + c->loci[i].key_off, key, len) == 0) {
      *locus = i; return FBS_GEN_OK;
    }
  return FBS_GEN_E_NOT_FOUND;
}

uint32_t fbs_gen_trait_count(const fbs_gen *c) { return ctx_ok(c) ? c->T : 0u; }
uint32_t fbs_gen_locus_count(const fbs_gen *c) { return ctx_ok(c) ? c->L : 0u; }

fbs_gen_status fbs_gen_set_season(fbs_gen *c, uint32_t season) {
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  if (season < c->season) return FBS_GEN_E_STATE;
  c->season = season;
  return FBS_GEN_OK;
}

uint32_t fbs_gen_season(const fbs_gen *c) { return ctx_ok(c) ? c->season : 0u; }
uint32_t fbs_gen_living(const fbs_gen *c) { return ctx_ok(c) ? c->living : 0u; }
uint32_t fbs_gen_record_count(const fbs_gen *c) { return ctx_ok(c) ? c->records : 0u; }

/* ------------------------------------------------------------------------- */
/* Records, living slots, phenotype                                          */
/* ------------------------------------------------------------------------- */

static uint32_t rec_find(const fbs_gen *c, uint32_t serial) {
  uint32_t lo = 0, hi = c->records;
  if (serial == 0u) return GEN_NONE;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2u;
    uint32_t s = c->rec[mid].serial;
    if (s == serial) return mid;
    if (s < serial) lo = mid + 1u; else hi = mid;
  }
  return GEN_NONE;
}

/* Living slot of id, or a status. */
static fbs_gen_status living_slot(const fbs_gen *c, fbs_gen_id id, uint32_t *slot) {
  uint32_t r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  if (!(c->rec[r].flags & FBS_GEN_ALIVE)) return FBS_GEN_E_DEAD;
  *slot = c->rec[r].slot;
  return FBS_GEN_OK;
}

static void phenotype(const fbs_gen *c, const int8_t *al, const int32_t *e, uint32_t f,
                      int32_t *z, int32_t *a) {
  int64_t num[FBS_GEN_MAX_TRAITS], add[FBS_GEN_MAX_TRAITS];
  uint32_t l, k, t;
  for (t = 0; t < c->T; t++) { num[t] = 0; add[t] = 0; }
  for (l = 0; l < c->L; l++) {
    const gen_locus *loc = &c->loci[l];
    int32_t a1 = al[l], a2 = al[c->L + l];
    int32_t hi = a1 > a2 ? a1 : a2, lo = a1 > a2 ? a2 : a1;
    for (k = loc->eff_begin; k < loc->eff_end; k++) {
      const gen_effect *ef = &c->effects[k];
      int64_t eh = effect_value(ef, hi), el = effect_value(ef, lo);
      num[ef->trait] += (int64_t)loc->h * eh + (int64_t)(256 - loc->h) * el;
      add[ef->trait] += eh + el;
    }
  }
  for (t = 0; t < c->T; t++) {
    const gen_trait *tr = &c->traits[t];
    int64_t v = (int64_t)tr->base + round_div(num[t], 128) + (int64_t)e[t]
              - round_div((int64_t)tr->depression * (int64_t)f, (int64_t)FBS_GEN_Q30);
    z[t] = clamp32(v, tr->lo, tr->hi);
    a[t] = sat32((int64_t)tr->base + add[t]);
  }
}

static void living_refresh(fbs_gen *c, uint32_t slot) {
  int32_t e[FBS_GEN_MAX_TRAITS], z[FBS_GEN_MAX_TRAITS], a[FBS_GEN_MAX_TRAITS];
  uint32_t t, cap = c->cfg.max_living;
  memset(e, 0, sizeof e);
  for (t = 0; t < c->T; t++) e[t] = c->l_e[(size_t)t * cap + slot];
  phenotype(c, c->l_al + (size_t)slot * c->stride, e, c->l_f[slot], z, a);
  for (t = 0; t < c->T; t++) { c->l_z[(size_t)t * cap + slot] = z[t]; c->l_a[(size_t)t * cap + slot] = a[t]; }
}

static void memo_clear(fbs_gen *c) { memset(c->memo, 0, sizeof(gen_memo) * c->cfg.memo_entries); }

/* Number of dead records eligible for eviction. */
static uint32_t evictable(const fbs_gen *c) {
  uint32_t i, n = 0;
  if (c->season < c->cfg.record_seasons) return 0;
  for (i = 0; i < c->records; i++)
    if (!(c->rec[i].flags & FBS_GEN_ALIVE) && c->rec[i].born < c->season - c->cfg.record_seasons) n++;
  return n;
}

/* Evict the k lowest-serial eligible dead records (k <= evictable). */
/* After compaction: parent record indices through rec_map (old -> new or GEN_NONE). */
static void remap_parents(fbs_gen *c) {
  uint32_t i;
  for (i = 0; i < c->records; i++) {
    gen_record *r = &c->rec[i];
    if (r->sire_r != GEN_NONE) r->sire_r = c->rec_map[r->sire_r];
    if (r->dam_r != GEN_NONE) r->dam_r = c->rec_map[r->dam_r];
  }
}

static void evict(fbs_gen *c, uint32_t k) {
  uint32_t i, w = 0, gone = 0, T = c->T;
  if (k == 0u) return;
  for (i = 0; i < c->records; i++) {
    const gen_record *r = &c->rec[i];
    if (gone < k && !(r->flags & FBS_GEN_ALIVE) && r->born < c->season - c->cfg.record_seasons) {
      gone++;
      c->rec_map[i] = GEN_NONE;
      continue;
    }
    c->rec_map[i] = w;
    if (w != i) {
      c->rec[w] = c->rec[i];
      if (c->cfg.record_traits) memmove(c->rec_z + (size_t)w * T, c->rec_z + (size_t)i * T, 4u * T);
    }
    w++;
  }
  c->records = w;
  remap_parents(c);
  memo_clear(c);
}

/* Capacity for `n` new individuals; evicts when needed. Returns evicted count via *ev. */
static fbs_gen_status reserve(fbs_gen *c, uint32_t n, uint32_t *ev, int apply) {
  uint32_t free_rec = c->cfg.max_records - c->records, need;
  if (n > c->cfg.max_living - c->living) return FBS_GEN_E_FULL;
  /* serials never wrap: next_serial 0 means every serial has been used */
  if (n > 0u && (c->next_serial == 0u || n - 1u > UINT32_MAX - c->next_serial)) return FBS_GEN_E_FULL;
  need = n > free_rec ? n - free_rec : 0u;
  if (need > 0u && evictable(c) < need) return FBS_GEN_E_FULL;
  if (apply) evict(c, need);
  if (ev) *ev = need;
  return FBS_GEN_OK;
}

/* Append a new individual: record + living slot. Caller reserved capacity. */
/* z/a: the phenotype and additive value already computed from (al, e, f), or NULL to compute. */
static uint32_t append_individual(fbs_gen *c, const int8_t *al, const int32_t *e, uint32_t f,
                                  uint32_t sire, uint32_t dam, uint32_t line, uint32_t flags,
                                  const int32_t *z, const int32_t *a) {
  uint32_t slot = c->living++, r = c->records++, t, cap = c->cfg.max_living;
  gen_record *rec = &c->rec[r];
  rec->serial = c->next_serial++;
  rec->sire = sire; rec->dam = dam; rec->born = c->season; rec->line = line; rec->f = f;
  rec->flags = FBS_GEN_ALIVE | flags; rec->children = 0; rec->honours = 0; rec->slot = slot;
  rec->sire_r = sire ? rec_find(c, sire) : GEN_NONE;
  rec->dam_r = dam ? rec_find(c, dam) : GEN_NONE;
  c->l_serial[slot] = rec->serial;
  c->l_f[slot] = f;
  c->l_stamp[slot] = 0;
  memset(c->l_al + (size_t)slot * c->stride, 0, c->stride);
  memcpy(c->l_al + (size_t)slot * c->stride, al, 2u * c->L);
  for (t = 0; t < c->T; t++) c->l_e[(size_t)t * cap + slot] = e[t];
  if (z && a) {
    for (t = 0; t < c->T; t++) { c->l_z[(size_t)t * cap + slot] = z[t]; c->l_a[(size_t)t * cap + slot] = a[t]; }
  } else {
    living_refresh(c, slot);
  }
  if (c->cfg.record_traits)
    for (t = 0; t < c->T; t++) c->rec_z[(size_t)r * c->T + t] = c->l_z[(size_t)t * cap + slot];
  return rec->serial;
}

static int64_t sd_eff(const fbs_gen *c, uint32_t t, uint32_t scale_q8) {
  return ((int64_t)c->traits[t].env_sd * (int64_t)scale_q8 + 128) >> 8;
}

fbs_gen_status fbs_gen_add_founders(fbs_gen *c, uint32_t count, uint32_t line, uint64_t seed,
                                    fbs_gen_id *out) {
  uint32_t i, l, s, t, j;
  int8_t al[2u * FBS_GEN_MAX_LOCI];
  int32_t e[FBS_GEN_MAX_TRAITS];
  fbs_gen_status st;
  gen_rng rng;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  if (count == 0u) return FBS_GEN_OK;
  st = reserve(c, count, NULL, 0);
  if (st != FBS_GEN_OK) return st;
  reserve(c, count, NULL, 1);
  for (i = 0; i < count; i++) {
    rng_stream(&rng, seed, i);
    for (s = 0; s < 2u; s++)
      for (l = 0; l < c->L; l++) {
        const gen_locus *loc = &c->loci[l];
        uint32_t r = rng_below(&rng, loc->fw_total);
        for (j = 0; j + 1u < loc->fcount && r >= loc->fw[j]; j++) r -= loc->fw[j];
        al[s * c->L + l] = loc->fa[j];
      }
    for (t = 0; t < c->T; t++) e[t] = rng_noise(&rng, sd_eff(c, t, 256u));
    {
      uint32_t id = append_individual(c, al, e, 0u, 0u, 0u, line, FBS_GEN_FOUNDER, NULL, NULL);
      if (out) out[i] = id;
    }
  }
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_add_genome(fbs_gen *c, const int8_t *alleles, const int32_t *nonheritable,
                                  uint32_t line, fbs_gen_id *out) {
  int32_t e[FBS_GEN_MAX_TRAITS];
  uint32_t l, t, id;
  fbs_gen_status st;
  if (!ctx_ok(c) || !alleles) return FBS_GEN_E_INVALID;
  for (l = 0; l < 2u * c->L; l++) {
    const gen_locus *loc = &c->loci[l % c->L];
    if (alleles[l] < loc->amin || alleles[l] > loc->amax) return FBS_GEN_E_RANGE;
  }
  for (t = 0; t < c->T; t++) {
    e[t] = nonheritable ? nonheritable[t] : 0;
    if (e[t] > (1 << 29) || e[t] < -(1 << 29)) return FBS_GEN_E_RANGE;
  }
  st = reserve(c, 1u, NULL, 0);
  if (st != FBS_GEN_OK) return st;
  reserve(c, 1u, NULL, 1);
  id = append_individual(c, alleles, e, 0u, 0u, 0u, line, FBS_GEN_FOUNDER, NULL, NULL);
  if (out) *out = id;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_retire(fbs_gen *c, fbs_gen_id id) {
  uint32_t r, slot, last, t, cap;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  if (!(c->rec[r].flags & FBS_GEN_ALIVE)) return FBS_GEN_E_DEAD;
  slot = c->rec[r].slot;
  last = c->living - 1u;
  cap = c->cfg.max_living;
  if (slot != last) {
    uint32_t mr = rec_find(c, c->l_serial[last]);
    c->l_serial[slot] = c->l_serial[last];
    c->l_f[slot] = c->l_f[last];
    c->l_stamp[slot] = c->l_stamp[last];
    memcpy(c->l_al + (size_t)slot * c->stride, c->l_al + (size_t)last * c->stride, c->stride);
    for (t = 0; t < c->T; t++) {
      size_t o = (size_t)t * cap;
      c->l_z[o + slot] = c->l_z[o + last];
      c->l_a[o + slot] = c->l_a[o + last];
      c->l_e[o + slot] = c->l_e[o + last];
    }
    c->rec[mr].slot = slot;
  }
  c->living--;
  c->rec[r].flags &= ~FBS_GEN_ALIVE;
  c->rec[r].slot = GEN_NONE;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_forget(fbs_gen *c, fbs_gen_id id) {
  uint32_t r, T;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  if (c->rec[r].flags & FBS_GEN_ALIVE) return FBS_GEN_E_STATE;
  T = c->T;
  {
    uint32_t i;
    for (i = 0; i < c->records; i++) c->rec_map[i] = i < r ? i : (i == r ? GEN_NONE : i - 1u);
  }
  memmove(&c->rec[r], &c->rec[r + 1u], sizeof(gen_record) * (c->records - r - 1u));
  if (c->cfg.record_traits)
    memmove(c->rec_z + (size_t)r * T, c->rec_z + (size_t)(r + 1u) * T, 4u * T * (c->records - r - 1u));
  c->records--;
  remap_parents(c);
  memo_clear(c);
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_set_line(fbs_gen *c, fbs_gen_id id, uint32_t line) {
  uint32_t r;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  c->rec[r].line = line;
  return FBS_GEN_OK;
}

static void info_of(const gen_record *r, fbs_gen_info *o) {
  o->id = r->serial; o->sire = r->sire; o->dam = r->dam; o->born = r->born; o->line = r->line;
  o->children = r->children; o->honours = r->honours; o->flags = r->flags; o->f_q30 = r->f;
}

fbs_gen_status fbs_gen_get(const fbs_gen *c, fbs_gen_id id, fbs_gen_info *out) {
  uint32_t r;
  if (!ctx_ok(c) || !out) return FBS_GEN_E_INVALID;
  r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  info_of(&c->rec[r], out);
  return FBS_GEN_OK;
}

static fbs_gen_status read_column_values(const fbs_gen *c, fbs_gen_id id, const int32_t *cols,
                                         int32_t *out, uint32_t cap) {
  uint32_t slot, t;
  fbs_gen_status st;
  if (!ctx_ok(c) || !out) return FBS_GEN_E_INVALID;
  st = living_slot(c, id, &slot);
  if (st != FBS_GEN_OK) return st;
  if (cap < c->T) return FBS_GEN_E_TRUNCATED;
  for (t = 0; t < c->T; t++) out[t] = cols[(size_t)t * c->cfg.max_living + slot];
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_traits(const fbs_gen *c, fbs_gen_id id, int32_t *z, uint32_t cap) {
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  return read_column_values(c, id, c->l_z, z, cap);
}
fbs_gen_status fbs_gen_additive(const fbs_gen *c, fbs_gen_id id, int32_t *a, uint32_t cap) {
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  return read_column_values(c, id, c->l_a, a, cap);
}
fbs_gen_status fbs_gen_nonheritable(const fbs_gen *c, fbs_gen_id id, int32_t *e, uint32_t cap) {
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  return read_column_values(c, id, c->l_e, e, cap);
}

fbs_gen_status fbs_gen_birth_traits(const fbs_gen *c, fbs_gen_id id, int32_t *z, uint32_t cap) {
  uint32_t r;
  if (!ctx_ok(c) || !z) return FBS_GEN_E_INVALID;
  if (!c->cfg.record_traits) return FBS_GEN_E_STATE;
  r = rec_find(c, id);
  if (r == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  if (cap < c->T) return FBS_GEN_E_TRUNCATED;
  memcpy(z, c->rec_z + (size_t)r * c->T, 4u * c->T);
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_genome(const fbs_gen *c, fbs_gen_id id, int8_t *alleles, uint32_t cap) {
  uint32_t slot;
  fbs_gen_status st;
  if (!ctx_ok(c) || !alleles) return FBS_GEN_E_INVALID;
  st = living_slot(c, id, &slot);
  if (st != FBS_GEN_OK) return st;
  if (cap < 2u * c->L) return FBS_GEN_E_TRUNCATED;
  memcpy(alleles, c->l_al + (size_t)slot * c->stride, 2u * c->L);
  return FBS_GEN_OK;
}

float fbs_gen_trait_f(const fbs_gen *c, fbs_gen_id id, uint32_t trait) {
  uint32_t slot;
  union { uint32_t u; float f; } nan_bits;
  nan_bits.u = 0x7FC00000u;
  if (!ctx_ok(c) || trait >= c->T || living_slot(c, id, &slot) != FBS_GEN_OK) return nan_bits.f;
  return (float)c->l_z[(size_t)trait * c->cfg.max_living + slot] / 1000.0f;
}

fbs_gen_status fbs_gen_column(const fbs_gen *c, uint32_t trait, const fbs_gen_id **ids,
                              const int32_t **values, uint32_t *count) {
  if (!ctx_ok(c) || !ids || !values || !count) return FBS_GEN_E_INVALID;
  if (trait >= c->T) return FBS_GEN_E_RANGE;
  *ids = c->l_serial;
  *values = c->l_z + (size_t)trait * c->cfg.max_living;
  *count = c->living;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_records(const fbs_gen *c, fbs_gen_info *out, size_t cap, size_t *count) {
  uint32_t i;
  if (!ctx_ok(c) || !count) return FBS_GEN_E_INVALID;
  *count = c->records;
  if (!out || cap == 0u) return c->records ? FBS_GEN_E_TRUNCATED : FBS_GEN_OK;
  if (cap < c->records) return FBS_GEN_E_TRUNCATED;
  for (i = 0; i < c->records; i++) info_of(&c->rec[i], &out[i]);
  return FBS_GEN_OK;
}

/* ------------------------------------------------------------------------- */
/* Kinship                                                                   */
/* ------------------------------------------------------------------------- */

static void sort_u32(uint32_t *v, uint32_t n) {
  /* heapsort ascending */
  uint32_t i, end;
  if (n < 2u) return;
  for (i = n / 2u; i-- > 0u;) {
    uint32_t root = i;
    for (;;) {
      uint32_t ch = 2u * root + 1u, tmp;
      if (ch >= n) break;
      if (ch + 1u < n && v[ch + 1u] > v[ch]) ch++;
      if (v[root] >= v[ch]) break;
      tmp = v[root]; v[root] = v[ch]; v[ch] = tmp; root = ch;
    }
  }
  for (end = n - 1u; end > 0u; end--) {
    uint32_t tmp = v[0], root = 0;
    v[0] = v[end]; v[end] = tmp;
    for (;;) {
      uint32_t ch = 2u * root + 1u;
      if (ch >= end) break;
      if (ch + 1u < end && v[ch + 1u] > v[ch]) ch++;
      if (v[root] >= v[ch]) break;
      tmp = v[root]; v[root] = v[ch]; v[ch] = tmp; root = ch;
    }
  }
}

/* Lays out the ancestor tree of id in heap order in k_*[side]; returns the count of known serials
 * copied (for the disjointness test) into list. */
static uint32_t kin_tree(fbs_gen *c, uint32_t side, uint32_t id, uint32_t *list) {
  uint32_t p, n = 0, H = c->k_heap, *ser = c->k_ser[side], *f = c->k_f[side], *ri = c->k_rec[side];
  uint8_t *has = c->k_has[side];
  ser[0] = id;
  ri[0] = rec_find(c, id);
  for (p = 0; p < H; p++) {
    uint32_t r = ser[p] ? ri[p] : GEN_NONE;
    has[p] = r != GEN_NONE;
    f[p] = has[p] ? c->rec[r].f : 0u;
    if (ser[p]) list[n++] = ser[p];
    if (2u * p + 2u < H) {
      ser[2u * p + 1u] = has[p] ? c->rec[r].sire : 0u;
      ser[2u * p + 2u] = has[p] ? c->rec[r].dam : 0u;
      ri[2u * p + 1u] = has[p] ? c->rec[r].sire_r : GEN_NONE;
      ri[2u * p + 2u] = has[p] ? c->rec[r].dam_r : GEN_NONE;
    }
  }
  return n;
}

static uint32_t heap_gen(uint32_t p) {
  uint32_t g = 0;
  p += 1u;
  while (p > 1u) { p >>= 1; g++; }
  return g;
}

/* Subtree masks of shared serials, bottom-up. */
static void kin_masks(fbs_gen *c, uint32_t side, uint32_t nshared) {
  uint32_t p = c->k_heap, H = c->k_heap;
  const uint32_t *ser = c->k_ser[side];
  uint64_t *mask = c->k_mask[side];
  while (p-- > 0u) {
    uint64_t m = 0;
    if (nshared > 64u) m = ~(uint64_t)0;
    else if (ser[p]) {
      uint32_t lo = 0, hi = nshared;
      while (lo < hi) {
        uint32_t mid = lo + (hi - lo) / 2u;
        if (c->k_shared[mid] < ser[p]) lo = mid + 1u; else hi = mid;
      }
      if (lo < nshared && c->k_shared[lo] == ser[p]) m = (uint64_t)1 << lo;
      if (2u * p + 2u < H) m |= mask[2u * p + 1u] | mask[2u * p + 2u];
    }
    mask[p] = m;
  }
}

/* f(x, y) for heap position x in tree 0 and y in tree 1: the recursion of the file comment. */
static uint32_t kin_dp(fbs_gen *c, uint32_t x, uint32_t y) {
  uint32_t sa = c->k_ser[0][x], sb = c->k_ser[1][y], v;
  uint64_t *m, sum;
  if (sa == 0u || sb == 0u) return 0u;
  if (sa == sb) return c->k_has[0][x] ? (uint32_t)(((uint64_t)FBS_GEN_Q30 + c->k_f[0][x]) >> 1) : 0u;
  if ((c->k_mask[0][x] & c->k_mask[1][y]) == 0u) return 0u;
  m = &c->k_memo[(size_t)x * c->k_heap + y];
  if ((uint32_t)(*m >> 32) == c->k_stamp) return (uint32_t)*m;
  if (sa > sb) {  /* x is younger: expand x */
    if (heap_gen(x) >= c->cfg.pedigree_depth || !c->k_has[0][x]) v = 0u;
    else {
      sum = (uint64_t)kin_dp(c, 2u * x + 1u, y) + (uint64_t)kin_dp(c, 2u * x + 2u, y);
      v = (uint32_t)(sum >> 1);
    }
  } else {        /* y is younger: expand y */
    if (heap_gen(y) >= c->cfg.pedigree_depth || !c->k_has[1][y]) v = 0u;
    else {
      sum = (uint64_t)kin_dp(c, x, 2u * y + 1u) + (uint64_t)kin_dp(c, x, 2u * y + 2u);
      v = (uint32_t)(sum >> 1);
    }
  }
  *m = ((uint64_t)c->k_stamp << 32) | v;
  return v;
}

static uint32_t memo_hash(uint32_t a, uint32_t b) {
  uint64_t k = ((uint64_t)a << 32) ^ ((uint64_t)b * 0x9E3779B97F4A7C15ULL);
  k ^= k >> 29; k *= 0xBF58476D1CE4E5B9ULL; k ^= k >> 32;
  return (uint32_t)k;
}

static uint32_t kinship(fbs_gen *c, uint32_t a, uint32_t b) {
  uint32_t na, nb, i = 0, j = 0, h, mask, v = 0, shared = 0;
  if (a == 0u || b == 0u) return 0u;
  if (a == b) {
    uint32_t r = rec_find(c, a);
    return r == GEN_NONE ? 0u : (uint32_t)(((uint64_t)FBS_GEN_Q30 + c->rec[r].f) >> 1);
  }
  if (a > b) { uint32_t t = a; a = b; b = t; }
  mask = c->cfg.memo_entries - 1u;
  h = memo_hash(a, b) & mask;
  for (i = 0; i < 8u; i++) {
    const gen_memo *m = &c->memo[(h + i) & mask];
    if (m->a == 0u) break;
    if (m->a == a && m->b == b) return m->v;
  }
  na = kin_tree(c, 0u, a, c->anc_a);
  nb = kin_tree(c, 1u, b, c->anc_b);
  sort_u32(c->anc_a, na);
  sort_u32(c->anc_b, nb);
  i = 0;
  while (i < na && j < nb) {
    if (c->anc_a[i] == c->anc_b[j]) {
      uint32_t x = c->anc_a[i];
      if (shared == 0u || (shared <= 64u && c->k_shared[shared - 1u] != x)) {
        if (shared < 64u) c->k_shared[shared] = x;
        shared++;
      }
      i++; j++;
    } else if (c->anc_a[i] < c->anc_b[j]) i++; else j++;
  }
  if (shared) {
    kin_masks(c, 0u, shared);
    kin_masks(c, 1u, shared);
    if (++c->k_stamp == 0u) {
      memset(c->k_memo, 0, 8u * (size_t)c->k_heap * c->k_heap);
      c->k_stamp = 1u;
    }
    v = kin_dp(c, 0u, 0u);
  }
  for (i = 0; i < 8u; i++) {
    gen_memo *m = &c->memo[(h + i) & mask];
    if (m->a == 0u) { m->a = a; m->b = b; m->g = 0u; m->v = v; break; }
  }
  return v;
}

fbs_gen_status fbs_gen_kinship(fbs_gen *c, fbs_gen_id a, fbs_gen_id b, uint32_t *q30) {
  if (!ctx_ok(c) || !q30) return FBS_GEN_E_INVALID;
  if (rec_find(c, a) == GEN_NONE || rec_find(c, b) == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  *q30 = kinship(c, a, b);
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_ancestors(const fbs_gen *c, fbs_gen_id id, uint32_t depth,
                                 fbs_gen_id *out, size_t cap, size_t *count) {
  size_t n, k;
  if (!ctx_ok(c) || !count) return FBS_GEN_E_INVALID;
  if (depth < 1u || depth > FBS_GEN_MAX_DEPTH) return FBS_GEN_E_RANGE;
  if (rec_find(c, id) == GEN_NONE) return FBS_GEN_E_NOT_FOUND;
  n = ((size_t)1 << (depth + 1u)) - 2u;
  *count = n;
  if (!out || cap < n) return FBS_GEN_E_TRUNCATED;
  for (k = 0; k < n; k++) {
    size_t pos = k + 2u, parent = pos / 2u;
    uint32_t child = parent == 1u ? id : out[parent - 2u], r;
    out[k] = 0u;
    if (child == 0u) continue;
    r = rec_find(c, child);
    if (r == GEN_NONE) continue;
    out[k] = (pos & 1u) ? c->rec[r].dam : c->rec[r].sire;
  }
  return FBS_GEN_OK;
}

/* ------------------------------------------------------------------------- */
/* Environment, indices, ranking, stats                                      */
/* ------------------------------------------------------------------------- */

void fbs_gen_environment_default(fbs_gen_environment *env) {
  if (!env) return;
  memset(env, 0, sizeof *env);
  env->mutation_scale_q8 = 256u;
  env->noise_scale_q8 = 256u;
}

static fbs_gen_status env_check(const fbs_gen *c, const fbs_gen_environment *env) {
  uint32_t i;
  if (!env) return FBS_GEN_OK;
  if (env->mutation_scale_q8 > FBS_GEN_SCALE_MAX_Q8 || env->noise_scale_q8 > FBS_GEN_SCALE_MAX_Q8)
    return FBS_GEN_E_RANGE;
  if (env->viability_count > FBS_GEN_MAX_VIABILITY) return FBS_GEN_E_RANGE;
  for (i = 0; i < FBS_GEN_MAX_TRAITS; i++)
    if (env->shift[i] > (1 << 28) || env->shift[i] < -(1 << 28)) return FBS_GEN_E_RANGE;
  for (i = 0; i < env->viability_count; i++) {
    const fbs_gen_viability *v = &env->viability[i];
    if (v->trait >= c->T || v->low > v->high) return FBS_GEN_E_RANGE;
  }
  return FBS_GEN_OK;
}

static fbs_gen_status index_check(const fbs_gen *c, const fbs_gen_index *ix) {
  (void)c;
  if (!ix) return FBS_GEN_E_INVALID;
  if (ix->basis > FBS_GEN_ADDITIVE || ix->standardize > 1u) return FBS_GEN_E_RANGE;
  return FBS_GEN_OK;
}

static int32_t living_value(const fbs_gen *c, uint32_t basis, uint32_t t, uint32_t slot) {
  const int32_t *col = basis == FBS_GEN_ADDITIVE ? c->l_a : c->l_z;
  return col[(size_t)t * c->cfg.max_living + slot];
}

static int64_t raw_score(const fbs_gen *c, const fbs_gen_index *ix, uint32_t slot) {
  int64_t s = 0;
  uint32_t t;
  for (t = 0; t < c->T; t++)
    if (ix->weight[t]) s += round_div((int64_t)ix->weight[t] * living_value(c, ix->basis, t, slot), 1000);
  return s;
}

/* Scores for n living slots. */
static void index_scores(const fbs_gen *c, const fbs_gen_index *ix, const uint32_t *slots, uint32_t n,
                         int64_t *out) {
  uint32_t i, t;
  for (i = 0; i < n; i++) out[i] = ix->standardize ? 0 : raw_score(c, ix, slots[i]);
  if (!ix->standardize || n == 0u) return;
  for (t = 0; t < c->T; t++) {
    int64_t sum = 0, mean;
    double var = 0.0, sd;
    int32_t sdi;
    if (!ix->weight[t]) continue;
    for (i = 0; i < n; i++) sum += living_value(c, ix->basis, t, slots[i]);
    mean = round_div(sum, (int64_t)n);
    for (i = 0; i < n; i++) {
      double d = (double)((int64_t)living_value(c, ix->basis, t, slots[i]) - mean);
      var += d * d;
    }
    sdi = sd_from_var(var / (double)n);
    sd = sdi > 0 ? (double)sdi : 1.0;
    for (i = 0; i < n; i++) {
      double d = (double)((int64_t)living_value(c, ix->basis, t, slots[i]) - mean);
      out[i] += round_d((double)ix->weight[t] * d / sd);
    }
  }
}

static int scored_before(const fbs_gen_scored *x, const fbs_gen_scored *y) {
  return x->score > y->score || (x->score == y->score && x->id < y->id);
}

static void sort_scored(fbs_gen_scored *v, uint32_t n) {
  /* heapsort so that "before" elements come first */
  uint32_t i, end;
  fbs_gen_scored tmp;
  if (n < 2u) return;
  for (i = n / 2u; i-- > 0u;) {
    uint32_t root = i;
    for (;;) {
      uint32_t ch = 2u * root + 1u;
      if (ch >= n) break;
      if (ch + 1u < n && scored_before(&v[ch], &v[ch + 1u])) ch++;
      if (!scored_before(&v[root], &v[ch])) break;
      tmp = v[root]; v[root] = v[ch]; v[ch] = tmp; root = ch;
    }
  }
  for (end = n - 1u; end > 0u; end--) {
    uint32_t root = 0;
    tmp = v[0]; v[0] = v[end]; v[end] = tmp;
    for (;;) {
      uint32_t ch = 2u * root + 1u;
      if (ch >= end) break;
      if (ch + 1u < end && scored_before(&v[ch], &v[ch + 1u])) ch++;
      if (!scored_before(&v[root], &v[ch])) break;
      tmp = v[root]; v[root] = v[ch]; v[ch] = tmp; root = ch;
    }
  }
}

/* Resolve ids to living slots (NULL/0 = all). */
static fbs_gen_status resolve_slots(const fbs_gen *c, const fbs_gen_id *ids, uint32_t n,
                                    uint32_t *slots, uint32_t *count) {
  uint32_t i;
  fbs_gen_status st;
  if (!ids || n == 0u) {
    for (i = 0; i < c->living; i++) slots[i] = i;
    *count = c->living;
    return FBS_GEN_OK;
  }
  if (n > c->cfg.max_living) return FBS_GEN_E_RANGE;
  for (i = 0; i < n; i++) {
    st = living_slot(c, ids[i], &slots[i]);
    if (st != FBS_GEN_OK) return st;
  }
  *count = n;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_rank(fbs_gen *c, const fbs_gen_index *index, const fbs_gen_id *candidates,
                            uint32_t n, fbs_gen_scored *out, size_t cap, size_t *count) {
  uint32_t m, i;
  fbs_gen_status st;
  if (!ctx_ok(c) || !count) return FBS_GEN_E_INVALID;
  st = index_check(c, index);
  if (st != FBS_GEN_OK) return st;
  st = resolve_slots(c, candidates, n, c->r_slot, &m);
  if (st != FBS_GEN_OK) return st;
  index_scores(c, index, c->r_slot, m, c->r_score);
  for (i = 0; i < m; i++) {
    c->r_buf[i].id = c->l_serial[c->r_slot[i]]; c->r_buf[i].pad = 0; c->r_buf[i].score = c->r_score[i];
  }
  sort_scored(c->r_buf, m);
  *count = m;
  if (out) for (i = 0; i < m && i < cap; i++) out[i] = c->r_buf[i];
  return FBS_GEN_OK;
}

/* Generic stats over rows: value of row k for trait t = z[k * stride]. */
typedef struct gen_rows {
  const int32_t *z, *a;
  size_t stride;
  const uint32_t *f;
  const uint32_t *idx;   /* row list or NULL for 0..n-1 */
  uint32_t n;
} gen_rows;

static void stats_rows(const gen_rows *rows, fbs_gen_trait_stats *o) {
  uint32_t k;
  int64_t sz = 0, sa = 0;
  uint64_t sf = 0;
  double vz = 0.0, va = 0.0, mz, ma;
  int32_t mn = INT32_MAX, mx = INT32_MIN;
  memset(o, 0, sizeof *o);
  o->n = rows->n;
  if (rows->n == 0u) return;
  for (k = 0; k < rows->n; k++) {
    size_t r = rows->idx ? rows->idx[k] : k;
    int32_t z = rows->z[r * rows->stride], a = rows->a[r * rows->stride];
    sz += z; sa += a; sf += rows->f[r];
    if (z < mn) mn = z;
    if (z > mx) mx = z;
  }
  mz = (double)sz / (double)rows->n;
  ma = (double)sa / (double)rows->n;
  for (k = 0; k < rows->n; k++) {
    size_t r = rows->idx ? rows->idx[k] : k;
    double dz = (double)rows->z[r * rows->stride] - mz, da = (double)rows->a[r * rows->stride] - ma;
    vz += dz * dz; va += da * da;
  }
  vz /= (double)rows->n; va /= (double)rows->n;
  o->mean = sat32(round_div(sz, (int64_t)rows->n));
  o->additive_mean = sat32(round_div(sa, (int64_t)rows->n));
  o->sd = sd_from_var(vz);
  o->additive_sd = sd_from_var(va);
  o->min = mn; o->max = mx;
  o->mean_f_q30 = (uint32_t)(sf / rows->n);
  if (vz > 0.0) {
    double h = va / vz;
    o->h2_q16 = h >= 1.0 ? 65536u : (uint32_t)round_d(h * 65536.0);
  }
}

fbs_gen_status fbs_gen_stats(const fbs_gen *c, uint32_t trait, const fbs_gen_id *subset,
                             uint32_t n, fbs_gen_trait_stats *out) {
  /* const: resolve slots on the fly in two passes via a small inline loop */
  gen_rows rows;
  uint32_t i, slot;
  fbs_gen_status st;
  int64_t sz = 0, sa = 0;
  uint64_t sf = 0;
  double vz = 0.0, va = 0.0, mz, ma;
  int32_t mn = INT32_MAX, mx = INT32_MIN;
  size_t cap;
  if (!ctx_ok(c) || !out) return FBS_GEN_E_INVALID;
  if (trait >= c->T) return FBS_GEN_E_RANGE;
  cap = c->cfg.max_living;
  if (!subset || n == 0u) {
    rows.z = c->l_z + (size_t)trait * cap; rows.a = c->l_a + (size_t)trait * cap; rows.stride = 1;
    rows.f = c->l_f; rows.idx = NULL; rows.n = c->living;
    stats_rows(&rows, out);
    return FBS_GEN_OK;
  }
  for (i = 0; i < n; i++) {
    st = living_slot(c, subset[i], &slot);
    if (st != FBS_GEN_OK) return st;
  }
  for (i = 0; i < n; i++) {
    int32_t z, a;
    living_slot(c, subset[i], &slot);
    z = c->l_z[(size_t)trait * cap + slot]; a = c->l_a[(size_t)trait * cap + slot];
    sz += z; sa += a; sf += c->l_f[slot];
    if (z < mn) mn = z;
    if (z > mx) mx = z;
  }
  mz = (double)sz / (double)n; ma = (double)sa / (double)n;
  for (i = 0; i < n; i++) {
    double dz, da;
    living_slot(c, subset[i], &slot);
    dz = (double)c->l_z[(size_t)trait * cap + slot] - mz;
    da = (double)c->l_a[(size_t)trait * cap + slot] - ma;
    vz += dz * dz; va += da * da;
  }
  vz /= (double)n; va /= (double)n;
  memset(out, 0, sizeof *out);
  out->n = n;
  out->mean = sat32(round_div(sz, (int64_t)n));
  out->additive_mean = sat32(round_div(sa, (int64_t)n));
  out->sd = sd_from_var(vz);
  out->additive_sd = sd_from_var(va);
  out->min = mn; out->max = mx;
  out->mean_f_q30 = (uint32_t)(sf / n);
  if (vz > 0.0) {
    double h = va / vz;
    out->h2_q16 = h >= 1.0 ? 65536u : (uint32_t)round_d(h * 65536.0);
  }
  return FBS_GEN_OK;
}

static fbs_gen_status check_subset(const fbs_gen *c, const fbs_gen_id *subset, uint32_t n, uint32_t *m) {
  uint32_t i, slot;
  fbs_gen_status st;
  if (!subset || n == 0u) { *m = c->living; return FBS_GEN_OK; }
  if (n > c->cfg.max_living) return FBS_GEN_E_RANGE;   /* keeps the Q16 sums in uint64 */
  for (i = 0; i < n; i++) {
    st = living_slot(c, subset[i], &slot);
    if (st != FBS_GEN_OK) return st;
  }
  *m = n;
  return FBS_GEN_OK;
}

static const int8_t *subset_genome(const fbs_gen *c, const fbs_gen_id *subset, uint32_t n, uint32_t k) {
  uint32_t slot = k;
  if (subset && n) living_slot(c, subset[k], &slot);
  return c->l_al + (size_t)slot * c->stride;
}

fbs_gen_status fbs_gen_diversity_of(const fbs_gen *c, const fbs_gen_id *subset, uint32_t n,
                                    fbs_gen_diversity *out) {
  uint32_t m, l, k, poly = 0;
  uint64_t het_sum = 0, exp_sum = 0;
  uint32_t counts[256];
  fbs_gen_status st;
  if (!ctx_ok(c) || !out) return FBS_GEN_E_INVALID;
  st = check_subset(c, subset, n, &m);
  if (st != FBS_GEN_OK) return st;
  memset(out, 0, sizeof *out);
  if (m == 0u) return FBS_GEN_OK;
  for (l = 0; l < c->L; l++) {
    uint64_t het = 0, sq = 0, total = 2ull * m;
    uint32_t distinct = 0, v;
    memset(counts, 0, sizeof counts);
    for (k = 0; k < m; k++) {
      const int8_t *g = subset_genome(c, subset, n, k);
      counts[(uint8_t)g[l]]++;
      counts[(uint8_t)g[c->L + l]]++;
      if (g[l] != g[c->L + l]) het++;
    }
    for (v = 0; v < 256u; v++)
      if (counts[v]) { distinct++; sq += (uint64_t)counts[v] * counts[v]; }
    if (distinct > 1u) poly++;
    het_sum += (het << 16) / m;
    exp_sum += 65536u - (sq << 16) / (total * total);
  }
  out->loci_polymorphic = poly;
  out->heterozygosity_q16 = (uint32_t)(het_sum / c->L);
  out->expected_heterozygosity_q16 = (uint32_t)(exp_sum / c->L);
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_allele_frequency(const fbs_gen *c, uint32_t locus, int8_t allele,
                                        const fbs_gen_id *subset, uint32_t n, uint32_t *q16) {
  uint32_t m, k;
  uint64_t hits = 0;
  fbs_gen_status st;
  if (!ctx_ok(c) || !q16) return FBS_GEN_E_INVALID;
  if (locus >= c->L) return FBS_GEN_E_RANGE;
  st = check_subset(c, subset, n, &m);
  if (st != FBS_GEN_OK) return st;
  if (m == 0u) { *q16 = 0; return FBS_GEN_OK; }
  for (k = 0; k < m; k++) {
    const int8_t *g = subset_genome(c, subset, n, k);
    hits += (uint64_t)(g[locus] == allele) + (uint64_t)(g[c->L + locus] == allele);
  }
  *q16 = (uint32_t)((hits << 16) / (2ull * m));
  return FBS_GEN_OK;
}

/* ------------------------------------------------------------------------- */
/* Breeding                                                                  */
/* ------------------------------------------------------------------------- */

static uint32_t pick_weighted(const uint64_t *cum, uint32_t n, uint64_t r) {
  uint32_t lo = 0, hi = n - 1u;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2u;
    if (cum[mid] > r) hi = mid; else lo = mid + 1u;
  }
  return lo;
}

fbs_gen_status fbs_gen_plan_pool(fbs_gen *c, const fbs_gen_pool_desc *d,
                                 const fbs_gen_contributor *con, uint32_t n,
                                 fbs_gen_union *out, uint32_t cap, uint32_t *written) {
  uint32_t i, k, w = 0;
  uint64_t total = 0;
  uint32_t kin_max;
  fbs_gen_status st;
  gen_rng rng;
  if (!ctx_ok(c) || !d || !con || !out || !written) return FBS_GEN_E_INVALID;
  if (d->children < 1u || d->children > c->cfg.max_brood) return FBS_GEN_E_RANGE;
  if (cap < d->children) return FBS_GEN_E_RANGE;
  if (n < 1u || n > c->cfg.max_contributors) return FBS_GEN_E_RANGE;
  if (d->redraws > 16u || d->kin_max_q16 > 65536u || d->allow_self > 1u) return FBS_GEN_E_RANGE;
  if (d->assort_k > 0u) {
    if (d->assort_k >= n && n > 1u) return FBS_GEN_E_RANGE;
    st = index_check(c, d->assort_index);
    if (st != FBS_GEN_OK) return st;
  }
  for (i = 0; i < n; i++) {
    st = living_slot(c, con[i].id, &c->p_slot[i]);
    if (st != FBS_GEN_OK) return st;
    total += con[i].weight;
    c->p_cum[i] = total;
  }
  if (total == 0u) return FBS_GEN_E_RANGE;
  kin_max = d->kin_max_q16 >= 65536u ? FBS_GEN_Q30 : d->kin_max_q16 << 14;
  if (d->assort_k > 0u) {
    index_scores(c, d->assort_index, c->p_slot, n, c->p_score);
    /* ascending score, ties by contributor index: sort_scored orders descending score, ascending id */
    for (i = 0; i < n; i++) { c->r_buf[i].id = i; c->r_buf[i].pad = 0; c->r_buf[i].score = -c->p_score[i]; }
    sort_scored(c->r_buf, n);
    for (i = 0; i < n; i++) { c->p_order[i] = c->r_buf[i].id; c->p_pos[c->r_buf[i].id] = i; }
  }
  rng_seed(&rng, d->seed);
  for (k = 0; k < d->children; k++) {
    uint32_t first = pick_weighted(c->p_cum, n, rng_below64(&rng, total));
    uint32_t best = GEN_NONE, best_kin = 0, t;
    for (t = 0; t <= d->redraws; t++) {
      uint32_t partner, kin;
      if (d->assort_k > 0u) {
        uint32_t p = c->p_pos[first], lo = p, hi = p + 1u, m = 0, j;
        int64_t s0 = c->p_score[first];
        uint64_t sum = 0, r;
        while (m < d->assort_k && (lo > 0u || hi < n)) {
          int take_lo;
          if (lo == 0u) take_lo = 0;
          else if (hi >= n) take_lo = 1;
          else {
            int64_t dl = s0 - c->p_score[c->p_order[lo - 1u]], dh = c->p_score[c->p_order[hi]] - s0;
            take_lo = dl <= dh;
          }
          if (take_lo) c->p_nb[m++] = c->p_order[--lo];
          else c->p_nb[m++] = c->p_order[hi++];
        }
        for (j = 0; j < m; j++) sum += con[c->p_nb[j]].weight;
        if (sum == 0u) continue;
        r = rng_below64(&rng, sum);
        for (j = 0; j + 1u < m && r >= con[c->p_nb[j]].weight; j++) r -= con[c->p_nb[j]].weight;
        partner = c->p_nb[j];
      } else {
        partner = pick_weighted(c->p_cum, n, rng_below64(&rng, total));
      }
      if (con[partner].id == con[first].id && !d->allow_self) continue;
      kin = kinship(c, con[first].id, con[partner].id);
      if (best == GEN_NONE || kin < best_kin) { best = partner; best_kin = kin; }
      if (kin <= kin_max) break;
    }
    if (best == GEN_NONE) continue;
    out[w].a = con[first].id;
    out[w].b = con[best].id;
    out[w].stream = d->stream_base + k;
    out[w].kinship_q30 = best_kin;
    w++;
  }
  *written = w;
  return FBS_GEN_OK;
}

static void make_gamete(const fbs_gen *c, const int8_t *parent, gen_rng *rng, int8_t *out,
                        uint32_t mscale, uint32_t *mutations) {
  uint32_t l, strand = rng_next(rng) >> 31;
  for (l = 0; l < c->L; l++) {
    const gen_locus *loc = &c->loci[l];
    int32_t a;
    uint32_t p;
    if (l > 0u && (rng_next(rng) >> 16) < loc->rec) strand ^= 1u;
    a = parent[strand * c->L + l];
    p = (uint32_t)(((uint64_t)loc->mut * mscale) >> 8);
    if (p > 65536u) p = 65536u;
    if (p > 0u && (rng_next(rng) >> 16) < p) {
      int32_t m = 1 + (int32_t)rng_below(rng, loc->step), old = a;
      if (rng_next(rng) >> 31) a += m; else a -= m;
      if (a < loc->amin) a = loc->amin;
      if (a > loc->amax) a = loc->amax;
      if (a != old) (*mutations)++;
    }
    out[l] = (int8_t)a;
  }
}

/* Validate and stage a brood in the scratch arrays. */
static fbs_gen_status brood_stage(fbs_gen *c, const fbs_gen_brood_desc *d, const fbs_gen_union *u,
                                  uint32_t n, fbs_gen_brood_report *rep, uint32_t *survivors) {
  fbs_gen_environment neutral;
  const fbs_gen_environment *env;
  uint32_t i, t, mscale, mutations = 0, alive = 0, T = c->T;
  int64_t sde[FBS_GEN_MAX_TRAITS];
  fbs_gen_status st;
  if (rep) rep->failed_union = n;
  if (!d || (n > 0u && !u)) return FBS_GEN_E_INVALID;
  if (n > c->cfg.max_brood) return FBS_GEN_E_RANGE;
  if (d->line_policy > FBS_GEN_LINE_HIGHER || d->allow_self > 1u) return FBS_GEN_E_RANGE;
  if (d->line_policy == FBS_GEN_LINE_HIGHER) {
    st = index_check(c, d->higher_index);
    if (st != FBS_GEN_OK) return st;
  }
  st = env_check(c, d->env);
  if (st != FBS_GEN_OK) return st;
  fbs_gen_environment_default(&neutral);
  env = d->env ? d->env : &neutral;
  for (i = 0; i < n; i++) {
    uint32_t sa, sb;
    st = living_slot(c, u[i].a, &sa);
    if (st == FBS_GEN_OK) st = living_slot(c, u[i].b, &sb);
    if (st == FBS_GEN_OK && u[i].a == u[i].b && !d->allow_self) st = FBS_GEN_E_INVALID;
    if (st != FBS_GEN_OK) { if (rep) rep->failed_union = i; return st; }
  }
  mscale = env->mutation_scale_q8;
  for (t = 0; t < T; t++) sde[t] = sd_eff(c, t, env->noise_scale_q8);
  for (i = 0; i < n; i++) {
    uint32_t sa = 0, sb = 0, v, f;
    int8_t *al = c->b_al + (size_t)i * c->stride;
    int32_t *e = c->b_e + (size_t)i * T;
    int32_t *z = c->b_z + (size_t)i * T;
    int32_t *a = c->b_a + (size_t)i * T;
    uint8_t ok = 1;
    gen_rng rng;
    living_slot(c, u[i].a, &sa);
    living_slot(c, u[i].b, &sb);
    f = kinship(c, u[i].a, u[i].b);
    rng_stream(&rng, d->seed, u[i].stream);
    make_gamete(c, c->l_al + (size_t)sa * c->stride, &rng, al, mscale, &mutations);
    make_gamete(c, c->l_al + (size_t)sb * c->stride, &rng, al + c->L, mscale, &mutations);
    for (t = 0; t < T; t++) e[t] = env->shift[t] + rng_noise(&rng, sde[t]);
    phenotype(c, al, e, f, z, a);
    for (v = 0; v < env->viability_count; v++) {
      const fbs_gen_viability *vr = &env->viability[v];
      int32_t x = z[vr->trait];
      uint32_t surv;
      if (x <= vr->low) surv = vr->floor_q16;
      else if (x >= vr->high) surv = 65536u;
      else surv = vr->floor_q16 + (uint32_t)(((uint64_t)(65536u - vr->floor_q16) * (uint64_t)((int64_t)x - vr->low))
                                             / (uint64_t)((int64_t)vr->high - vr->low));
      if ((rng_next(&rng) >> 16) >= surv) ok = 0;
    }
    c->b_f[i] = f;
    c->b_ok[i] = ok;
    alive += ok;
    switch (d->line_policy) {
      case FBS_GEN_LINE_FIRST: c->b_line[i] = c->rec[rec_find(c, u[i].a)].line; break;
      case FBS_GEN_LINE_HIGHER: {
        int64_t xa = raw_score(c, d->higher_index, sa), xb = raw_score(c, d->higher_index, sb);
        c->b_line[i] = c->rec[rec_find(c, xa >= xb ? u[i].a : u[i].b)].line;
        break;
      }
      default: c->b_line[i] = d->line; break;
    }
  }
  if (rep) {
    uint64_t sf = 0;
    for (i = 0; i < n; i++) if (c->b_ok[i]) sf += c->b_f[i];
    rep->born = alive; rep->culled = n - alive; rep->mutations = mutations; rep->evicted_records = 0;
    rep->mean_f_q30 = alive ? (uint32_t)(sf / alive) : 0u;
    rep->failed_union = n;
  }
  *survivors = alive;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_brood(fbs_gen *c, const fbs_gen_brood_desc *d, const fbs_gen_union *u,
                             uint32_t n, fbs_gen_id *out_children, fbs_gen_brood_report *report) {
  fbs_gen_brood_report rep;
  uint32_t survivors = 0, ev = 0, i;
  fbs_gen_status st;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  memset(&rep, 0, sizeof rep);
  st = brood_stage(c, d, u, n, &rep, &survivors);
  if (st != FBS_GEN_OK) { if (report) report->failed_union = rep.failed_union; return st; }
  st = reserve(c, survivors, &ev, 0);
  if (st != FBS_GEN_OK) { if (report) report->failed_union = n; return st; }
  reserve(c, survivors, &ev, 1);
  rep.evicted_records = ev;
  c->brood_counter++;
  for (i = 0; i < n; i++) {
    uint32_t id = 0, p;
    if (c->b_ok[i])
      id = append_individual(c, c->b_al + (size_t)i * c->stride, c->b_e + (size_t)i * c->T, c->b_f[i],
                             u[i].a, u[i].b, c->b_line[i], 0u,
                             c->b_z + (size_t)i * c->T, c->b_a + (size_t)i * c->T);
    for (p = 0; p < 2u; p++) {
      fbs_gen_id pid = p ? u[i].b : u[i].a;
      uint32_t r = rec_find(c, pid), slot;
      if (p == 1u && pid == u[i].a) break;
      slot = c->rec[r].slot;
      if (c->l_stamp[slot] != c->brood_counter) { c->l_stamp[slot] = c->brood_counter; c->rec[r].honours++; }
      if (id) c->rec[r].children++;
    }
    if (out_children) out_children[i] = id;
  }
  if (report) *report = rep;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_brood_preview(fbs_gen *c, const fbs_gen_brood_desc *d, const fbs_gen_union *u,
                                     uint32_t n, fbs_gen_trait_stats *out, uint32_t cap,
                                     fbs_gen_brood_report *report) {
  fbs_gen_brood_report rep;
  uint32_t survivors = 0, i, k = 0, t;
  fbs_gen_status st;
  if (!ctx_ok(c)) return FBS_GEN_E_INVALID;
  memset(&rep, 0, sizeof rep);
  st = brood_stage(c, d, u, n, &rep, &survivors);
  if (st != FBS_GEN_OK) { if (report) report->failed_union = rep.failed_union; return st; }
  if (report) *report = rep;
  if (!out) return FBS_GEN_OK;
  /* survivor row list; b_line is not needed after staging */
  for (i = 0; i < n; i++) if (c->b_ok[i]) c->b_line[k++] = i;
  for (t = 0; t < c->T && t < cap; t++) {
    gen_rows rows;
    rows.z = c->b_z + t; rows.a = c->b_a + t; rows.stride = c->T;
    rows.f = c->b_f; rows.idx = c->b_line; rows.n = k;
    stats_rows(&rows, &out[t]);
  }
  return cap < c->T ? FBS_GEN_E_TRUNCATED : FBS_GEN_OK;
}

/* Raw forecast of one union: mean (before clamping) and variance per trait. */
static void forecast_raw(fbs_gen *c, uint32_t sa, uint32_t sb, uint32_t f,
                         const fbs_gen_environment *env, int64_t *mean, double *var) {
  int64_t num[FBS_GEN_MAX_TRAITS];
  int64_t tv[FBS_GEN_MAX_TRAITS][4];
  uint8_t touched[FBS_GEN_MAX_TRAITS];
  const int8_t *A = c->l_al + (size_t)sa * c->stride, *B = c->l_al + (size_t)sb * c->stride;
  uint32_t l, t, k, e;
  for (t = 0; t < c->T; t++) { num[t] = 0; var[t] = 0.0; touched[t] = 0; tv[t][0] = tv[t][1] = tv[t][2] = tv[t][3] = 0; }
  for (l = 0; l < c->L; l++) {
    const gen_locus *loc = &c->loci[l];
    if (loc->eff_begin == loc->eff_end) continue;
    for (k = 0; k < 4u; k++) {
      int32_t x = A[(k & 1u) * c->L + l], y = B[(k >> 1) * c->L + l];
      int32_t hi = x > y ? x : y, lo = x > y ? y : x;
      for (e = loc->eff_begin; e < loc->eff_end; e++) {
        const gen_effect *ef = &c->effects[e];
        tv[ef->trait][k] += (int64_t)loc->h * effect_value(ef, hi) + (int64_t)(256 - loc->h) * effect_value(ef, lo);
        touched[ef->trait] = 1;
      }
    }
    for (t = 0; t < c->T; t++) {
      double s, ss;
      if (!touched[t]) continue;
      s = (double)(tv[t][0] + tv[t][1] + tv[t][2] + tv[t][3]);
      ss = (double)tv[t][0] * (double)tv[t][0] + (double)tv[t][1] * (double)tv[t][1]
         + (double)tv[t][2] * (double)tv[t][2] + (double)tv[t][3] * (double)tv[t][3];
      num[t] += tv[t][0] + tv[t][1] + tv[t][2] + tv[t][3];
      var[t] += (ss / 4.0 - (s / 4.0) * (s / 4.0)) / 16384.0;
      touched[t] = 0; tv[t][0] = tv[t][1] = tv[t][2] = tv[t][3] = 0;
    }
  }
  for (t = 0; t < c->T; t++) {
    const gen_trait *tr = &c->traits[t];
    double sd = (double)sd_eff(c, t, env->noise_scale_q8);
    mean[t] = (int64_t)tr->base + round_div(num[t], 512) + env->shift[t]
            - round_div((int64_t)tr->depression * (int64_t)f, (int64_t)FBS_GEN_Q30);
    var[t] += sd * sd;
  }
}

static void forecast_out(const fbs_gen *c, uint32_t t, int64_t mean, double var, fbs_gen_forecast *o) {
  const gen_trait *tr = &c->traits[t];
  int32_t sd = sd_from_var(var);
  int64_t half = round_div((int64_t)sd * 83991, 65536); /* 1.2816 in Q16 */
  o->mean = clamp32(mean, tr->lo, tr->hi);
  o->sd = sd;
  o->p10 = clamp32(mean - half, tr->lo, tr->hi);
  o->p90 = clamp32(mean + half, tr->lo, tr->hi);
}

fbs_gen_status fbs_gen_forecast_union(fbs_gen *c, fbs_gen_id a, fbs_gen_id b,
                                      const fbs_gen_environment *env,
                                      fbs_gen_forecast *out, uint32_t cap, uint32_t *f_q30) {
  fbs_gen_environment neutral;
  int64_t mean[FBS_GEN_MAX_TRAITS];
  double var[FBS_GEN_MAX_TRAITS];
  uint32_t sa = 0u, sb = 0u, f, t; /* set by living_slot; zeroed so MSVC's C4701 flow check agrees */
  fbs_gen_status st;
  if (!ctx_ok(c) || !out) return FBS_GEN_E_INVALID;
  st = env_check(c, env);
  if (st != FBS_GEN_OK) return st;
  st = living_slot(c, a, &sa);
  if (st == FBS_GEN_OK) st = living_slot(c, b, &sb);
  if (st != FBS_GEN_OK) return st;
  if (cap < c->T) return FBS_GEN_E_TRUNCATED;
  fbs_gen_environment_default(&neutral);
  f = kinship(c, a, b);
  forecast_raw(c, sa, sb, f, env ? env : &neutral, mean, var);
  for (t = 0; t < c->T; t++) forecast_out(c, t, mean[t], var[t], &out[t]);
  if (f_q30) *f_q30 = f;
  return FBS_GEN_OK;
}

fbs_gen_status fbs_gen_forecast_plan(fbs_gen *c, const fbs_gen_union *u, uint32_t n,
                                     const fbs_gen_environment *env,
                                     fbs_gen_forecast *out, uint32_t cap) {
  fbs_gen_environment neutral;
  int64_t mean[FBS_GEN_MAX_TRAITS];
  double var[FBS_GEN_MAX_TRAITS], sm[FBS_GEN_MAX_TRAITS], smm[FBS_GEN_MAX_TRAITS], sv[FBS_GEN_MAX_TRAITS];
  int64_t msum[FBS_GEN_MAX_TRAITS];
  uint32_t i, t, sa, sb;
  fbs_gen_status st;
  if (!ctx_ok(c) || !out || !u) return FBS_GEN_E_INVALID;
  if (n < 1u) return FBS_GEN_E_RANGE;
  st = env_check(c, env);
  if (st != FBS_GEN_OK) return st;
  for (i = 0; i < n; i++) {
    st = living_slot(c, u[i].a, &sa);
    if (st == FBS_GEN_OK) st = living_slot(c, u[i].b, &sb);
    if (st != FBS_GEN_OK) return st;
  }
  if (cap < c->T) return FBS_GEN_E_TRUNCATED;
  fbs_gen_environment_default(&neutral);
  for (t = 0; t < c->T; t++) { sm[t] = 0.0; smm[t] = 0.0; sv[t] = 0.0; msum[t] = 0; }
  for (i = 0; i < n; i++) {
    living_slot(c, u[i].a, &sa);
    living_slot(c, u[i].b, &sb);
    forecast_raw(c, sa, sb, kinship(c, u[i].a, u[i].b), env ? env : &neutral, mean, var);
    for (t = 0; t < c->T; t++) {
      msum[t] += mean[t];
      sm[t] += (double)mean[t];
      smm[t] += (double)mean[t] * (double)mean[t];
      sv[t] += var[t];
    }
  }
  for (t = 0; t < c->T; t++) {
    double m = sm[t] / (double)n, v = sv[t] / (double)n + (smm[t] / (double)n - m * m);
    if (v < 0.0) v = 0.0;
    forecast_out(c, t, round_div(msum[t], (int64_t)n), v, &out[t]);
  }
  return FBS_GEN_OK;
}

/* ------------------------------------------------------------------------- */
/* Snapshots                                                                 */
/* ------------------------------------------------------------------------- */

static size_t living_bytes(const fbs_gen *c) { return 4u + 4u * (size_t)c->T + 2u * (size_t)c->L; }
static size_t record_bytes(const fbs_gen *c) { return GEN_RECORD_BYTES + (c->cfg.record_traits ? 4u * (size_t)c->T : 0u); }

size_t fbs_gen_snapshot_size(const fbs_gen *c) {
  if (!ctx_ok(c)) return 0;
  return GEN_HEADER_BYTES + (size_t)c->living * living_bytes(c) + (size_t)c->records * record_bytes(c) + 8u;
}

fbs_gen_status fbs_gen_save(const fbs_gen *c, void *out, size_t cap, size_t *bytes) {
  size_t need;
  uint8_t *p;
  uint32_t i, t;
  if (!ctx_ok(c) || !bytes) return FBS_GEN_E_INVALID;
  need = fbs_gen_snapshot_size(c);
  *bytes = need;
  if (!out || cap < need) return FBS_GEN_E_TRUNCATED;
  p = (uint8_t *)out;
  memcpy(p, "FBSGEN\0\1", 8);
  put_u32(p + 8, FBS_GEN_SNAPSHOT_VERSION);
  put_u32(p + 12, FBS_GEN_VERSION);
  put_u64(p + 16, c->hash_structure);
  put_u64(p + 24, c->hash_tuning);
  put_u32(p + 32, c->L);
  put_u32(p + 36, c->T);
  put_u32(p + 40, c->cfg.record_traits);
  put_u32(p + 44, c->season);
  put_u32(p + 48, c->next_serial);
  put_u32(p + 52, c->living);
  put_u32(p + 56, c->records);
  put_u32(p + 60, 0u);
  p += GEN_HEADER_BYTES;
  for (i = 0; i < c->living; i++) {
    put_u32(p, c->l_serial[i]); p += 4;
    for (t = 0; t < c->T; t++) { put_u32(p, (uint32_t)c->l_e[(size_t)t * c->cfg.max_living + i]); p += 4; }
    memcpy(p, c->l_al + (size_t)i * c->stride, 2u * c->L); p += 2u * c->L;
  }
  for (i = 0; i < c->records; i++) {
    const gen_record *r = &c->rec[i];
    put_u32(p, r->serial); put_u32(p + 4, r->sire); put_u32(p + 8, r->dam); put_u32(p + 12, r->born);
    put_u32(p + 16, r->line); put_u32(p + 20, r->f); put_u32(p + 24, r->flags);
    put_u32(p + 28, r->children); put_u32(p + 32, r->honours);
    p += GEN_RECORD_BYTES;
    if (c->cfg.record_traits)
      for (t = 0; t < c->T; t++) { put_u32(p, (uint32_t)c->rec_z[(size_t)i * c->T + t]); p += 4; }
  }
  put_u64(p, fnv_bytes(FNV_INIT, out, need - 8u));
  return FBS_GEN_OK;
}

/* Binary search a record serial in the image. */
static uint32_t image_rec_find(const uint8_t *recs, size_t rb, uint32_t count, uint32_t serial) {
  uint32_t lo = 0, hi = count;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2u, s = get_u32(recs + (size_t)mid * rb);
    if (s == serial) return mid;
    if (s < serial) lo = mid + 1u; else hi = mid;
  }
  return GEN_NONE;
}

fbs_gen_status fbs_gen_load(fbs_gen *c, const void *data, size_t bytes, uint32_t flags) {
  const uint8_t *p = (const uint8_t *)data, *liv, *recs;
  uint32_t L, T, rt, season, next, living, records, i, t, alive = 0, prev = 0;
  size_t lb, rb, need;
  if (!ctx_ok(c) || !data) return FBS_GEN_E_INVALID;
  if (flags & ~FBS_GEN_LOAD_RETUNE) return FBS_GEN_E_RANGE;
  if (bytes < GEN_HEADER_BYTES + 8u) return FBS_GEN_E_CORRUPT;
  if (memcmp(p, "FBSGEN\0\1", 8) != 0) return FBS_GEN_E_CORRUPT;
  if (get_u64(p + bytes - 8u) != fnv_bytes(FNV_INIT, p, bytes - 8u)) return FBS_GEN_E_CORRUPT;
  if (get_u32(p + 8) != FBS_GEN_SNAPSHOT_VERSION) return FBS_GEN_E_SCHEMA;
  L = get_u32(p + 32); T = get_u32(p + 36); rt = get_u32(p + 40);
  if (L != c->L || T != c->T || get_u64(p + 16) != c->hash_structure) return FBS_GEN_E_SCHEMA;
  if (get_u64(p + 24) != c->hash_tuning && !(flags & FBS_GEN_LOAD_RETUNE)) return FBS_GEN_E_SCHEMA;
  if (rt != c->cfg.record_traits) return FBS_GEN_E_SCHEMA;
  season = get_u32(p + 44); next = get_u32(p + 48); living = get_u32(p + 52); records = get_u32(p + 56);
  if (get_u32(p + 60) != 0u) return FBS_GEN_E_CORRUPT;
  if (living > c->cfg.max_living || records > c->cfg.max_records) return FBS_GEN_E_FULL;
  if (living > records) return FBS_GEN_E_CORRUPT;   /* next 0: every serial used */
  lb = living_bytes(c); rb = record_bytes(c);
  need = GEN_HEADER_BYTES + (size_t)living * lb + (size_t)records * rb + 8u;
  if (bytes != need) return FBS_GEN_E_CORRUPT;
  liv = p + GEN_HEADER_BYTES;
  recs = liv + (size_t)living * lb;
  /* pass one: records */
  for (i = 0; i < records; i++) {
    const uint8_t *r = recs + (size_t)i * rb;
    uint32_t serial = get_u32(r), sire = get_u32(r + 4), dam = get_u32(r + 8), born = get_u32(r + 12);
    uint32_t f = get_u32(r + 20), fl = get_u32(r + 24);
    if (serial == 0u || serial <= prev || (next != 0u && serial >= next)) return FBS_GEN_E_CORRUPT;
    if (sire >= serial || dam >= serial || born > season || f > FBS_GEN_Q30) return FBS_GEN_E_CORRUPT;
    if (fl & ~(FBS_GEN_ALIVE | FBS_GEN_FOUNDER)) return FBS_GEN_E_CORRUPT;
    if (fl & FBS_GEN_ALIVE) alive++;
    prev = serial;
  }
  if (alive != living) return FBS_GEN_E_CORRUPT;
  /* living: each maps to an alive record, no duplicates, alleles in range */
  for (i = 0; i < living; i++) {
    const uint8_t *e = liv + (size_t)i * lb;
    uint32_t serial = get_u32(e), r = image_rec_find(recs, rb, records, serial), l;
    if (r == GEN_NONE || !(get_u32(recs + (size_t)r * rb + 24) & FBS_GEN_ALIVE)) return FBS_GEN_E_CORRUPT;
    for (t = 0; t < T; t++) {
      int32_t v = (int32_t)get_u32(e + 4 + 4u * t);
      if (v > (1 << 29) || v < -(1 << 29)) return FBS_GEN_E_CORRUPT;
    }
    for (l = 0; l < 2u * L; l++) {
      int8_t a = (int8_t)e[4u + 4u * T + l];
      const gen_locus *loc = &c->loci[l % L];
      if (a < loc->amin || a > loc->amax) return FBS_GEN_E_CORRUPT;
    }
    c->r_slot[i] = serial;   /* scratch only */
  }
  sort_u32(c->r_slot, living);
  for (i = 1; i < living; i++) if (c->r_slot[i] == c->r_slot[i - 1u]) return FBS_GEN_E_CORRUPT;
  /* pass two: decode */
  c->season = season; c->next_serial = next; c->records = records; c->living = living;
  c->brood_counter = 0;
  for (i = 0; i < records; i++) {
    const uint8_t *r = recs + (size_t)i * rb;
    gen_record *d = &c->rec[i];
    d->serial = get_u32(r); d->sire = get_u32(r + 4); d->dam = get_u32(r + 8); d->born = get_u32(r + 12);
    d->line = get_u32(r + 16); d->f = get_u32(r + 20); d->flags = get_u32(r + 24);
    d->children = get_u32(r + 28); d->honours = get_u32(r + 32); d->slot = GEN_NONE;
    if (rt) for (t = 0; t < T; t++) c->rec_z[(size_t)i * T + t] = (int32_t)get_u32(r + GEN_RECORD_BYTES + 4u * t);
  }
  for (i = 0; i < records; i++) {
    gen_record *d = &c->rec[i];
    d->sire_r = d->sire ? rec_find(c, d->sire) : GEN_NONE;
    d->dam_r = d->dam ? rec_find(c, d->dam) : GEN_NONE;
  }
  for (i = 0; i < living; i++) {
    const uint8_t *e = liv + (size_t)i * lb;
    uint32_t r = rec_find(c, get_u32(e));
    c->l_serial[i] = get_u32(e);
    c->l_f[i] = c->rec[r].f;
    c->l_stamp[i] = 0;
    c->rec[r].slot = i;
    for (t = 0; t < T; t++) c->l_e[(size_t)t * c->cfg.max_living + i] = (int32_t)get_u32(e + 4 + 4u * t);
    memset(c->l_al + (size_t)i * c->stride, 0, c->stride);
    memcpy(c->l_al + (size_t)i * c->stride, e + 4u + 4u * T, 2u * L);
    living_refresh(c, i);
  }
  memo_clear(c);
  return FBS_GEN_OK;
}
