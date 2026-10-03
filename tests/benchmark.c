/* tests/benchmark.c — performance spot check for fbs_genetics on the goblin preset (design.md
 * section 12): 20,000 living, 32,768 records, pedigree depth 5, a population bred for several
 * seasons so kinship has real pedigrees to walk. Prints the median of repeated runs per operation.
 * Host code (POSIX clock_gettime). Build with -DFBS_GENETICS_BUILD_BENCHMARK=ON, Release, idle machine.
 */
#define _POSIX_C_SOURCE 199309L
#include "fbs/genetics.h"
#include "goblin.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define REPS 21

static double now_us(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1e3;
}

static int cmp_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x > y;
}

static double median(double *v, int n) {
  qsort(v, (size_t)n, sizeof *v, cmp_double);
  return v[n / 2];
}

int main(void) {
  goblin_schema_data *sd = (goblin_schema_data *)malloc(sizeof *sd);
  fbs_gen_config cfg = fbs_gen_config_default();
  fbs_gen *g = NULL;
  void *mem;
  size_t bytes, img_n, got, cnt;
  uint8_t *img, *img2;
  fbs_gen_id *ids = (fbs_gen_id *)malloc(sizeof *ids * 20000u), *kids = (fbs_gen_id *)malloc(sizeof *kids * 2048u);
  fbs_gen_scored *top = (fbs_gen_scored *)malloc(sizeof *top * 20000u);
  fbs_gen_contributor con[400];
  fbs_gen_union *un = (fbs_gen_union *)malloc(sizeof *un * 2048u);
  fbs_gen_pool_desc pd;
  fbs_gen_brood_desc bd;
  fbs_gen_brood_report rep;
  fbs_gen_index ix;
  fbs_gen_trait_stats st[GOB_TRAITS];
  fbs_gen_forecast fc[GOB_TRAITS];
  double t[REPS], t0;
  uint32_t i, s, written, living, q, related_a = 0, related_b = 0;
  int r;
  const fbs_gen_id *col; const int32_t *vals;

  if (!sd || !ids || !kids || !top || !un) return 1;
  goblin_schema_build(sd, 40);
  cfg.max_living = 20000; cfg.max_records = 32768; cfg.pedigree_depth = 5; cfg.max_brood = 2048;
  cfg.max_contributors = 400; cfg.record_seasons = 8;
  bytes = fbs_gen_memory_for(&cfg, &sd->schema);
  mem = malloc(bytes);
  if (!mem || fbs_gen_create(&cfg, &sd->schema, mem, bytes, &g) != FBS_GEN_OK) return 1;
  printf("memory_for: %.2f MB (goblin schema, 20,000 living, 32,768 records)\n", (double)bytes / 1048576.0);

  memset(&ix, 0, sizeof ix);
  ix.weight[GOB_SIZE] = 1000; ix.weight[GOB_STRENGTH] = 1000; ix.standardize = 1;
  fbs_gen_add_founders(g, 8000, 0, 1, NULL);
  /* breed until the population is full and the pedigree is several generations deep */
  for (s = 1; s <= 10u; s++) {
    fbs_gen_set_season(g, s);
    while (fbs_gen_living(g) + 2000u > 20000u) {
      fbs_gen_column(g, GOB_SIZE, &col, &vals, &living);
      fbs_gen_retire(g, col[0]);
    }
    fbs_gen_rank(g, &ix, NULL, 0, top, 400, &cnt);
    for (i = 0; i < 400u; i++) { con[i].id = top[i].id; con[i].weight = 1000u; }
    memset(&pd, 0, sizeof pd);
    pd.children = 2000; pd.seed = s; pd.kin_max_q16 = 8192; pd.redraws = 4;
    fbs_gen_plan_pool(g, &pd, con, 400, un, 2048, &written);
    memset(&bd, 0, sizeof bd); bd.seed = 100u + s;
    { int bs = fbs_gen_brood(g, &bd, un, written, kids, &rep); if (bs != FBS_GEN_OK) { printf("setup brood failed: %s\n", fbs_gen_status_name(bs)); return 1; } }
    related_a = kids[0]; related_b = kids[1];
  }
  while (fbs_gen_living(g) + 2000u > 20000u) {
    fbs_gen_column(g, GOB_SIZE, &col, &vals, &living);
    fbs_gen_retire(g, col[0]);
  }
  fbs_gen_set_season(g, 11);
  printf("population: %u living, %u records, season 11\n", fbs_gen_living(g), fbs_gen_record_count(g));

  img_n = fbs_gen_snapshot_size(g);
  img = (uint8_t *)malloc(img_n); img2 = (uint8_t *)malloc(img_n);
  if (!img || !img2) return 1;
  fbs_gen_save(g, img, img_n, &got);

  fbs_gen_rank(g, &ix, NULL, 0, top, 400, &cnt);
  for (i = 0; i < 400u; i++) { con[i].id = top[i].id; con[i].weight = 1000u + i; }
  memset(&pd, 0, sizeof pd);
  pd.children = 2000; pd.seed = 77; pd.kin_max_q16 = 8192; pd.redraws = 4;

  for (r = 0; r < REPS; r++) {
    fbs_gen_load(g, img, img_n, 0);          /* cold memo each time */
    t0 = now_us(); fbs_gen_plan_pool(g, &pd, con, 400, un, 2048, &written); t[r] = now_us() - t0;
  }
  printf("plan_pool: 2,000 unions from 400 contributors, kin redraw 4 (cold memo): %8.1f us\n", median(t, REPS));
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_plan_pool(g, &pd, con, 400, un, 2048, &written); t[r] = now_us() - t0;
  }
  printf("plan_pool: same, warm memo:                                    %8.1f us\n", median(t, REPS));
  memset(&bd, 0, sizeof bd); bd.seed = 99;
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_brood_preview(g, &bd, un, written, st, GOB_TRAITS, &rep); t[r] = now_us() - t0;
  }
  printf("brood_preview: %u children:                                  %8.1f us\n", written, median(t, REPS));
  for (r = 0; r < REPS; r++) {
    fbs_gen_load(g, img, img_n, 0);
    t0 = now_us(); fbs_gen_brood(g, &bd, un, written, kids, &rep); t[r] = now_us() - t0;
  }
  printf("brood: %u children, cold memo (born %u, evicted %u records): %8.1f us\n", written, rep.born,
         rep.evicted_records, median(t, REPS));
  for (r = 0; r < REPS; r++) {
    fbs_gen_load(g, img, img_n, 0);
    fbs_gen_plan_pool(g, &pd, con, 400, un, 2048, &written);   /* the game's order: plan, then commit */
    t0 = now_us(); fbs_gen_brood(g, &bd, un, written, kids, &rep); t[r] = now_us() - t0;
  }
  printf("brood: %u children right after its plan (warm memo):        %8.1f us\n", written, median(t, REPS));
  fbs_gen_load(g, img, img_n, 0);
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_rank(g, &ix, NULL, 0, top, 20000, &cnt); t[r] = now_us() - t0;
  }
  printf("rank: %u living, 2 traits, standardized, full sort:        %8.1f us\n", (unsigned)cnt, median(t, REPS));
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_forecast_union(g, un[r].a, un[r].b, NULL, fc, GOB_TRAITS, &q); t[r] = now_us() - t0;
  }
  printf("forecast_union:                                              %8.2f us\n", median(t, REPS));
  for (r = 0; r < REPS; r++) {
    fbs_gen_load(g, img, img_n, 0);
    t0 = now_us(); fbs_gen_kinship(g, ids[0] ? ids[0] : 1u, related_a, &q); t[r] = now_us() - t0;
  }
  printf("kinship: founder vs newborn (fast path, cold):               %8.2f us\n", median(t, REPS));
  for (r = 0; r < REPS; r++) {
    fbs_gen_load(g, img, img_n, 0);
    t0 = now_us(); fbs_gen_kinship(g, related_a, related_b, &q); t[r] = now_us() - t0;
  }
  printf("kinship: two newborns of one pool (related, cold memo), f %.4f: %6.2f us\n", (double)q / FBS_GEN_Q30,
         median(t, REPS));
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_save(g, img2, img_n, &got); t[r] = now_us() - t0;
  }
  printf("save: %.2f MB image:                                          %8.1f us\n", (double)got / 1048576.0,
         median(t, REPS));
  for (r = 0; r < REPS; r++) {
    t0 = now_us(); fbs_gen_load(g, img, img_n, 0); t[r] = now_us() - t0;
  }
  printf("load (validate + decode + recompute phenotypes):              %8.1f us\n", median(t, REPS));
  (void)ids;
  free(img); free(img2); free(mem); free(sd); free(ids); free(kids); free(top); free(un);
  return 0;
}
