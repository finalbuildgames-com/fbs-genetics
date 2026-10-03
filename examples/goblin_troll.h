/* examples/goblin_troll.h — the prototype's troll-line gathering (design.md section 8, config A by
 * default) written against the public API over the goblin preset.
 *
 * Host code, not part of the library: it allocates with malloc. Shared by tests/test_genetics.c
 * and examples/gathering.c so the assertion and the printed table run the same simulation.
 */
#ifndef FBS_GENETICS_EXAMPLES_GOBLIN_TROLL_H
#define FBS_GENETICS_EXAMPLES_GOBLIN_TROLL_H

#include "goblin.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Troll line: design.md section 8, config A by default                      */
/* ------------------------------------------------------------------------- */

#define GOB_TROLL_LINE 1u
#define GOB_MAX_SEASONS 16

typedef struct goblin_run_config {
  uint32_t seasons;          /* <= GOB_MAX_SEASONS */
  uint32_t founders;         /* 1000 */
  uint32_t found_n;          /* 30 set aside at the first gathering */
  uint32_t honour_n;         /* 20 honoured in the troll pen */
  uint32_t brood;            /* 60 x mean fecundity */
  uint32_t kin_max_q16;      /* 8192 = 0.125 */
  uint32_t redraws;          /* 4 */
  uint32_t workers_honoured; /* 200 */
  uint32_t workers_brood;    /* 400 */
  uint32_t worker_cap;       /* 1000 */
  uint32_t lifespan;         /* 4 */
  uint16_t giant_permille;   /* 40 */
} goblin_run_config;

static goblin_run_config goblin_run_config_a(void) {
  goblin_run_config c;
  c.seasons = 12; c.founders = 1000; c.found_n = 30; c.honour_n = 20; c.brood = 60;
  c.kin_max_q16 = 8192; c.redraws = 4; c.workers_honoured = 200; c.workers_brood = 400;
  c.worker_cap = 1000; c.lifespan = 4; c.giant_permille = 40;
  return c;
}

typedef struct goblin_row {
  uint32_t season, line_n, brood;
  int32_t size, strength, deftness, stride, climb, stamina, fecundity; /* line means, milli */
  uint32_t f_q30;           /* line mean F */
  uint32_t troll_permille;  /* share of the line with size >= 2.5 */
  int32_t worker_deftness;
  uint32_t worker_f_q30;
  uint32_t giant1_q16, giant2_q16; /* line allele frequencies */
} goblin_row;

typedef struct goblin_run_result {
  goblin_row rows[GOB_MAX_SEASONS];
  uint32_t first_season_2500;  /* first season with line mean size >= 2.5, 0 = never */
  uint64_t snapshot_digest;    /* FNV-1a 64 of the final snapshot */
  size_t snapshot_bytes;
  int status;
} goblin_run_result;

static uint64_t goblin_fnv(const uint8_t *p, size_t n) {
  uint64_t h = 0xCBF29CE484222325ULL;
  size_t i;
  for (i = 0; i < n; i++) { h ^= p[i]; h *= 0x100000001B3ULL; }
  return h;
}

/* Runs the troll line. When keep is non-NULL the context survives the run: *keep receives it and
 * *keep_mem its block, which the caller frees. */
static int goblin_run_troll_keep(uint64_t seed, const goblin_run_config *rc, goblin_run_result *res,
                                 fbs_gen **keep, void **keep_mem) {
  goblin_schema_data *sd = (goblin_schema_data *)malloc(sizeof *sd);
  fbs_gen_config cfg = fbs_gen_config_default();
  fbs_gen *g = NULL;
  void *mem = NULL;
  size_t bytes, max_serial;
  uint8_t *age = NULL;
  fbs_gen_id *ids = NULL, *tc = NULL, *wc = NULL, *kids = NULL;
  fbs_gen_scored *top = NULL;
  fbs_gen_contributor *con = NULL;
  fbs_gen_union *un = NULL;
  fbs_gen_index troll_ix, worker_ix;
  uint32_t giant1 = 0, giant2 = 0, s, i;
  int st = FBS_GEN_OK;
  memset(res, 0, sizeof *res);
  if (!sd) return FBS_GEN_E_MEMORY;
  goblin_schema_build(sd, rc->giant_permille);
  cfg.max_living = 4096; cfg.max_records = 8192; cfg.max_brood = 1024; cfg.max_contributors = 256;
  bytes = fbs_gen_memory_for(&cfg, &sd->schema);
  max_serial = 1u + rc->founders + (size_t)rc->seasons * (2u * rc->brood + rc->workers_brood) + 16u;
  mem = malloc(bytes);
  age = (uint8_t *)calloc(max_serial, 1);
  ids = (fbs_gen_id *)malloc(sizeof *ids * cfg.max_living);
  tc = (fbs_gen_id *)malloc(sizeof *tc * cfg.max_living);
  wc = (fbs_gen_id *)malloc(sizeof *wc * cfg.max_living);
  kids = (fbs_gen_id *)malloc(sizeof *kids * cfg.max_brood);
  top = (fbs_gen_scored *)malloc(sizeof *top * cfg.max_living);
  con = (fbs_gen_contributor *)malloc(sizeof *con * cfg.max_contributors);
  un = (fbs_gen_union *)malloc(sizeof *un * cfg.max_brood);
  if (!bytes || !mem || !age || !ids || !tc || !wc || !kids || !top || !con || !un) { st = FBS_GEN_E_MEMORY; goto done; }
  st = fbs_gen_create(&cfg, &sd->schema, mem, bytes, &g);
  if (st) goto done;
  fbs_gen_locus_find(g, "GIANT1", 6, &giant1);
  fbs_gen_locus_find(g, "GIANT2", 6, &giant2);
  memset(&troll_ix, 0, sizeof troll_ix);
  troll_ix.weight[GOB_SIZE] = 1000; troll_ix.weight[GOB_STRENGTH] = 1000; troll_ix.standardize = 1;
  memset(&worker_ix, 0, sizeof worker_ix);
  worker_ix.weight[GOB_DEFTNESS] = 1000; worker_ix.weight[GOB_STAMINA] = 1000; worker_ix.standardize = 1;

  st = fbs_gen_add_founders(g, rc->founders, 0u, seed, ids);
  if (st) goto done;
  for (i = 0; i < rc->founders; i++) age[ids[i]] = 1;
  { /* first gathering: set aside the biggest and strongest as the troll line */
    size_t n = 0;
    st = fbs_gen_rank(g, &troll_ix, NULL, 0, top, rc->found_n, &n);
    if (st) goto done;
    for (i = 0; i < rc->found_n && i < n; i++) fbs_gen_set_line(g, top[i].id, GOB_TROLL_LINE);
  }

  for (s = 1; s <= rc->seasons; s++) {
    const fbs_gen_id *col_ids; const int32_t *col_size;
    uint32_t living, ntc = 0, nwc = 0, k, written;
    fbs_gen_brood_desc bd;
    fbs_gen_pool_desc pd;
    fbs_gen_brood_report rep;
    goblin_row *row = &res->rows[s - 1u];
    fbs_gen_set_season(g, s);
    fbs_gen_column(g, GOB_SIZE, &col_ids, &col_size, &living);
    memcpy(ids, col_ids, sizeof *ids * living);
    for (i = 0; i < living; i++) {
      fbs_gen_info info;
      uint32_t need = col_size[i] >= 2000 ? 2u : 1u;
      fbs_gen_get(g, ids[i], &info);
      if (age[ids[i]] < need) continue;
      if (info.line == GOB_TROLL_LINE) tc[ntc++] = ids[i]; else wc[nwc++] = ids[i];
    }
    memset(&bd, 0, sizeof bd);
    bd.line_policy = FBS_GEN_LINE_EXPLICIT;
    /* troll pen */
    row->brood = 0;
    if (ntc >= 2u) {
      size_t n = 0;
      int64_t fsum = 0;
      int32_t z[GOB_TRAITS];
      uint32_t honoured, nb;
      st = fbs_gen_rank(g, &troll_ix, tc, ntc, top, rc->honour_n, &n);
      if (st) goto done;
      honoured = (uint32_t)(n < rc->honour_n ? n : rc->honour_n);
      for (k = 0; k < honoured; k++) {
        fbs_gen_traits(g, top[k].id, z, GOB_TRAITS);
        con[k].id = top[k].id; con[k].weight = (uint32_t)z[GOB_FECUNDITY];
        fsum += z[GOB_FECUNDITY];
      }
      {
        int64_t mf = fsum / honoured;
        if (mf < 100) mf = 100;
        if (mf > 2000) mf = 2000;
        nb = (uint32_t)((rc->brood * (uint64_t)mf + 500u) / 1000u);
      }
      if (honoured >= 2u && nb > 0u) {
        memset(&pd, 0, sizeof pd);
        pd.children = nb; pd.seed = seed * 1000003u + s * 2u; pd.kin_max_q16 = rc->kin_max_q16;
        pd.redraws = rc->redraws;
        st = fbs_gen_plan_pool(g, &pd, con, honoured, un, cfg.max_brood, &written);
        if (st) goto done;
        bd.seed = seed * 7919u + s * 2u; bd.line = GOB_TROLL_LINE;
        st = fbs_gen_brood(g, &bd, un, written, kids, &rep);
        if (st) goto done;
        for (k = 0; k < written; k++) if (kids[k]) age[kids[k]] = 0;
        row->brood = rep.born;
      }
    }
    /* workers */
    if (nwc >= 2u) {
      size_t n = 0;
      int32_t z[GOB_TRAITS];
      uint32_t honoured;
      st = fbs_gen_rank(g, &worker_ix, wc, nwc, top, rc->workers_honoured, &n);
      if (st) goto done;
      honoured = (uint32_t)(n < rc->workers_honoured ? n : rc->workers_honoured);
      for (k = 0; k < honoured; k++) {
        fbs_gen_traits(g, top[k].id, z, GOB_TRAITS);
        con[k].id = top[k].id; con[k].weight = (uint32_t)z[GOB_FECUNDITY];
      }
      memset(&pd, 0, sizeof pd);
      pd.children = rc->workers_brood; pd.seed = seed * 1000003u + s * 2u + 1u; pd.kin_max_q16 = 4096;
      pd.redraws = 2;
      st = fbs_gen_plan_pool(g, &pd, con, honoured, un, cfg.max_brood, &written);
      if (st) goto done;
      bd.seed = seed * 7919u + s * 2u + 1u; bd.line = 0u;
      st = fbs_gen_brood(g, &bd, un, written, kids, &rep);
      if (st) goto done;
      for (k = 0; k < written; k++) if (kids[k]) age[kids[k]] = 0;
    }
    /* ageing and death */
    fbs_gen_column(g, GOB_SIZE, &col_ids, &col_size, &living);
    memcpy(ids, col_ids, sizeof *ids * living);
    for (i = 0; i < living; i++) {
      age[ids[i]]++;
      if (age[ids[i]] > rc->lifespan) fbs_gen_retire(g, ids[i]);
    }
    { /* worker cap: oldest die first, ties by serial */
      uint32_t nw = 0, a;
      fbs_gen_column(g, GOB_SIZE, &col_ids, &col_size, &living);
      memcpy(ids, col_ids, sizeof *ids * living);
      for (i = 0; i < living; i++) {
        fbs_gen_info info;
        fbs_gen_get(g, ids[i], &info);
        if (info.line != GOB_TROLL_LINE) wc[nw++] = ids[i];
      }
      if (nw > rc->worker_cap) {
        uint32_t excess = nw - rc->worker_cap;
        for (a = rc->lifespan + 1u; a-- > 0u && excess > 0u;) {
          /* serials ascending within an age: wc is not sorted, so scan by serial */
          uint32_t last = 0;
          for (;;) {
            uint32_t best = 0;
            for (i = 0; i < nw; i++)
              if (wc[i] > last && age[wc[i]] == a && (best == 0u || wc[i] < best)) best = wc[i];
            if (best == 0u || excess == 0u) break;
            fbs_gen_retire(g, best);
            excess--;
            last = best;
          }
        }
      }
    }
    { /* line row */
      fbs_gen_trait_stats ts;
      uint32_t nl = 0, nw = 0, trolls = 0;
      int32_t z[GOB_TRAITS];
      fbs_gen_column(g, GOB_SIZE, &col_ids, &col_size, &living);
      memcpy(ids, col_ids, sizeof *ids * living);
      for (i = 0; i < living; i++) {
        fbs_gen_info info;
        fbs_gen_get(g, ids[i], &info);
        if (info.line == GOB_TROLL_LINE) {
          tc[nl++] = ids[i];
          fbs_gen_traits(g, ids[i], z, GOB_TRAITS);
          if (z[GOB_SIZE] >= 2500) trolls++;
        } else {
          wc[nw++] = ids[i];
        }
      }
      row->season = s; row->line_n = nl;
      if (nl > 0u) {
        fbs_gen_stats(g, GOB_SIZE, tc, nl, &ts); row->size = ts.mean; row->f_q30 = ts.mean_f_q30;
        fbs_gen_stats(g, GOB_STRENGTH, tc, nl, &ts); row->strength = ts.mean;
        fbs_gen_stats(g, GOB_DEFTNESS, tc, nl, &ts); row->deftness = ts.mean;
        fbs_gen_stats(g, GOB_STRIDE, tc, nl, &ts); row->stride = ts.mean;
        fbs_gen_stats(g, GOB_CLIMB, tc, nl, &ts); row->climb = ts.mean;
        fbs_gen_stats(g, GOB_STAMINA, tc, nl, &ts); row->stamina = ts.mean;
        fbs_gen_stats(g, GOB_FECUNDITY, tc, nl, &ts); row->fecundity = ts.mean;
        row->troll_permille = trolls * 1000u / nl;
        fbs_gen_allele_frequency(g, giant1, 1, tc, nl, &row->giant1_q16);
        fbs_gen_allele_frequency(g, giant2, 1, tc, nl, &row->giant2_q16);
      }
      if (nw > 0u) {
        fbs_gen_stats(g, GOB_DEFTNESS, wc, nw, &ts);
        row->worker_deftness = ts.mean; row->worker_f_q30 = ts.mean_f_q30;
      }
      if (!res->first_season_2500 && row->size >= 2500) res->first_season_2500 = s;
    }
  }
  {
    size_t n = fbs_gen_snapshot_size(g), got = 0;
    uint8_t *img = (uint8_t *)malloc(n);
    if (!img) { st = FBS_GEN_E_MEMORY; goto done; }
    st = fbs_gen_save(g, img, n, &got);
    if (st == FBS_GEN_OK) { res->snapshot_digest = goblin_fnv(img, got); res->snapshot_bytes = got; }
    free(img);
  }
done:
  res->status = st;
  if (keep && keep_mem && st == FBS_GEN_OK) { *keep = g; *keep_mem = mem; mem = NULL; }
  free(un); free(con); free(top); free(kids); free(wc); free(tc); free(ids); free(age); free(mem); free(sd);
  return st;
}

#define goblin_run_troll(seed, rc, res) goblin_run_troll_keep((seed), (rc), (res), NULL, NULL)

#endif /* FBS_GENETICS_EXAMPLES_GOBLIN_TROLL_H */
