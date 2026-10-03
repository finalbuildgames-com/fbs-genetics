/* tests/test_genetics_stats.c — statistical tests for fbs_genetics, public header only, fixed seeds.
 *
 * What it proves:
 *   mendel      Aa x Aa at an unlinked locus gives 1:2:1 genotypes and, with h = 256, 3:1
 *               phenotypes (chi-square below the p = 0.001 critical value), 3 seeds x 40,000;
 *   linkage     the recombinant fraction between two loci is within 3 SE of r for r = 0.01, 0.1,
 *               0.3 and 0.5, 3 seeds x 20,000 gametes;
 *   breeder     one generation of 20% truncation selection on an additive unlinked schema: the
 *               realized response R is close to h^2 S per seed and on average (5 seeds);
 *   forecast    the analytic union forecast's mean is within 3 SE of 20,000 committed children and
 *               its SD within 5% (unlinked schema with dominance); the goblin schema's (linked)
 *               errors are measured and printed; forecast_plan of one union equals forecast_union;
 *   troll line  the goblin preset's config-A troll line over 12 seasons, 3 seeds, lands in the
 *               prototype's range (design.md section 8); the table is printed.
 */
#include "fbs/genetics.h"
#include "goblin_troll.h"

#include <math.h>
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
#define CHECK_OK(expr)                                                          \
  do {                                                                          \
    int status_ = (int)(expr);                                                  \
    CHECK(status_ == FBS_GEN_OK, "%s = %s", #expr, fbs_gen_status_name(status_)); \
  } while (0)

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

static void locus(fbs_gen_locus_desc *l, const char *key, int8_t amin, int8_t amax, uint16_t h, uint16_t rec) {
  memset(l, 0, sizeof *l);
  l->key = key; l->key_len = strlen(key);
  l->amin = amin; l->amax = amax; l->h = h; l->recombination_q16 = rec;
  l->founder_count = 1; l->founder[0].allele = amin; l->founder[0].weight = 1;
}

static fbs_gen_config big_config(uint32_t living, uint32_t brood) {
  fbs_gen_config cfg = fbs_gen_config_default();
  cfg.max_living = living; cfg.max_records = living; cfg.max_brood = brood; cfg.max_contributors = 64;
  return cfg;
}

/* ---- Mendelian segregation and dominance --------------------------------- */

static void test_mendel(void) {
  static const fbs_gen_trait_desc traits[2] = {
    {"add", 3, 0, -10000, 10000, 0, 0}, {"dom", 3, 0, -10000, 10000, 0, 0} };
  fbs_gen_locus_desc loci[2];
  fbs_gen_effect effects[2] = { {0, 0, 1, 1000}, {1, 1, 1, 1000} };
  fbs_gen_schema s;
  fbs_gen_config cfg = big_config(40004, 40000);
  fbs_gen_union *un = (fbs_gen_union *)malloc(sizeof *un * 40000u);
  fbs_gen_id *kids = (fbs_gen_id *)malloc(sizeof *kids * 40000u);
  uint32_t seed, i;
  locus(&loci[0], "A", 0, 1, 128, 32768);
  locus(&loci[1], "D", 0, 1, 256, 32768);
  s.traits = traits; s.trait_count = 2; s.loci = loci; s.locus_count = 2; s.effects = effects; s.effect_count = 2;
  if (!un || !kids) { CHECK(0, "alloc"); free(un); free(kids); return; }
  for (seed = 1; seed <= 3u; seed++) {
    void *mem;
    fbs_gen *g = make_ctx(&cfg, &s, &mem);
    int8_t al[4] = { 0, 0, 1, 1 };          /* a1 = (0, 0), a2 = (1, 1): Aa at both loci */
    fbs_gen_id p, q;
    fbs_gen_brood_desc bd;
    fbs_gen_brood_report rep;
    uint32_t geno[3] = { 0, 0, 0 }, dom = 0;
    double chi = 0.0, chi_d;
    if (!g) { CHECK(0, "mendel ctx"); continue; }
    CHECK_OK(fbs_gen_add_genome(g, al, NULL, 0, &p));
    CHECK_OK(fbs_gen_add_genome(g, al, NULL, 0, &q));
    for (i = 0; i < 40000u; i++) { un[i].a = p; un[i].b = q; un[i].stream = i; un[i].kinship_q30 = 0; }
    memset(&bd, 0, sizeof bd); bd.seed = seed * 1000003u;
    CHECK_OK(fbs_gen_brood(g, &bd, un, 40000, kids, &rep));
    for (i = 0; i < 40000u; i++) {
      int32_t z[2];
      int8_t gk[4];
      fbs_gen_traits(g, kids[i], z, 2);
      fbs_gen_genome(g, kids[i], gk, 4);
      geno[gk[0] + gk[2]]++;
      if (z[1] == 2000) dom++;
      CHECK(z[0] == 1000 * (gk[0] + gk[2]), "additive phenotype");
      CHECK(z[1] == ((gk[1] | gk[3]) ? 2000 : 0), "dominant phenotype");
    }
    chi = ((double)geno[0] - 10000.0) * ((double)geno[0] - 10000.0) / 10000.0
        + ((double)geno[1] - 20000.0) * ((double)geno[1] - 20000.0) / 20000.0
        + ((double)geno[2] - 10000.0) * ((double)geno[2] - 10000.0) / 10000.0;
    chi_d = ((double)dom - 30000.0) * ((double)dom - 30000.0) / 30000.0
          + ((double)(40000u - dom) - 10000.0) * ((double)(40000u - dom) - 10000.0) / 10000.0;
    printf("mendel seed %u: aa %u Aa %u AA %u chi2 %.2f; dominant %u recessive %u chi2 %.2f\n", seed,
           geno[0], geno[1], geno[2], chi, dom, 40000u - dom, chi_d);
    CHECK(chi < 13.82, "1:2:1 chi2 %.2f", chi);
    CHECK(chi_d < 10.83, "3:1 chi2 %.2f", chi_d);
    free(mem);
  }
  free(un); free(kids);
}

/* ---- linkage ------------------------------------------------------------- */

static void test_linkage(void) {
  static const fbs_gen_trait_desc trait = {"z", 1, 0, 0, 0, 0, 0};
  static const uint16_t r_q16[4] = { 655, 6554, 19661, 32768 };
  static const char *keys[8] = { "A0", "B0", "A1", "B1", "A2", "B2", "A3", "B3" };
  fbs_gen_locus_desc loci[8];
  fbs_gen_schema s;
  fbs_gen_config cfg = big_config(20004, 20000);
  fbs_gen_union *un = (fbs_gen_union *)malloc(sizeof *un * 20000u);
  fbs_gen_id *kids = (fbs_gen_id *)malloc(sizeof *kids * 20000u);
  uint32_t seed, i, k;
  for (k = 0; k < 4u; k++) {
    locus(&loci[2u * k], keys[2u * k], 0, 1, 128, 32768);
    locus(&loci[2u * k + 1u], keys[2u * k + 1u], 0, 1, 128, r_q16[k]);
  }
  s.traits = &trait; s.trait_count = 1; s.loci = loci; s.locus_count = 8; s.effects = NULL; s.effect_count = 0;
  if (!un || !kids) { CHECK(0, "alloc"); free(un); free(kids); return; }
  for (seed = 1; seed <= 3u; seed++) {
    void *mem;
    fbs_gen *g = make_ctx(&cfg, &s, &mem);
    int8_t cis[16], tester[16];
    fbs_gen_id p, q;
    fbs_gen_brood_desc bd;
    uint32_t rec[4] = { 0, 0, 0, 0 };
    if (!g) { CHECK(0, "linkage ctx"); continue; }
    memset(cis, 0, sizeof cis); memset(tester, 0, sizeof tester);
    for (i = 0; i < 8u; i++) cis[i] = 1;     /* AB/ab: strand 1 all ones, strand 2 all zeros */
    CHECK_OK(fbs_gen_add_genome(g, cis, NULL, 0, &p));
    CHECK_OK(fbs_gen_add_genome(g, tester, NULL, 0, &q));
    for (i = 0; i < 20000u; i++) { un[i].a = p; un[i].b = q; un[i].stream = i; un[i].kinship_q30 = 0; }
    memset(&bd, 0, sizeof bd); bd.seed = seed * 7u + 1u;
    CHECK_OK(fbs_gen_brood(g, &bd, un, 20000, kids, NULL));
    for (i = 0; i < 20000u; i++) {
      int8_t gk[16];
      fbs_gen_genome(g, kids[i], gk, 16);
      for (k = 0; k < 4u; k++) rec[k] += gk[2u * k] != gk[2u * k + 1u];
    }
    for (k = 0; k < 4u; k++) {
      double r = (double)r_q16[k] / 65536.0, obs = (double)rec[k] / 20000.0;
      double se = sqrt(r * (1.0 - r) / 20000.0);
      printf("linkage seed %u: r %.4f observed %.4f (%.2f SE)\n", seed, r, obs, (obs - r) / se);
      CHECK(fabs(obs - r) <= 3.0 * se, "r %.4f observed %.4f", r, obs);
    }
    free(mem);
  }
  free(un); free(kids);
}

/* ---- an additive unlinked schema for breeder and forecast tests ------------ */

typedef struct quant_schema {
  fbs_gen_trait_desc trait[2];
  fbs_gen_locus_desc loci[40];
  fbs_gen_effect effects[80];
  char keys[40][8];
  fbs_gen_schema schema;
} quant_schema;

/* `loci` unlinked loci, alleles -8..8 with founders uniform -2..2, weight w on trait 0; trait 1
 * reads the same loci with weight w/2 and dominance h. env_sd per trait. */
static void quant_build(quant_schema *q, uint32_t loci, int32_t w, int32_t env_sd, uint16_t h) {
  uint32_t i, j, ne = 0;
  memset(q, 0, sizeof *q);
  q->trait[0].key = "z"; q->trait[0].key_len = 1; q->trait[0].base = 10000;
  q->trait[0].lo = -1000000; q->trait[0].hi = 1000000; q->trait[0].env_sd = env_sd;
  q->trait[1] = q->trait[0]; q->trait[1].key = "y";
  for (i = 0; i < loci; i++) {
    fbs_gen_locus_desc *l = &q->loci[i];
    snprintf(q->keys[i], sizeof q->keys[i], "Q%u", (unsigned)i % 100u);
    l->key = q->keys[i]; l->key_len = strlen(q->keys[i]);
    l->amin = -8; l->amax = 8; l->h = i % 2u ? h : 128; l->recombination_q16 = 32768;
    l->founder_count = 5;
    for (j = 0; j < 5u; j++) { l->founder[j].allele = (int8_t)((int)j - 2); l->founder[j].weight = 1; }
    q->effects[ne].locus = (uint16_t)i; q->effects[ne].trait = 0; q->effects[ne].match = (int8_t)FBS_GEN_ANY;
    q->effects[ne].weight = w; ne++;
    q->effects[ne].locus = (uint16_t)i; q->effects[ne].trait = 1; q->effects[ne].match = (int8_t)FBS_GEN_ANY;
    q->effects[ne].weight = w / 2; ne++;
  }
  q->schema.traits = q->trait; q->schema.trait_count = 2;
  q->schema.loci = q->loci; q->schema.locus_count = loci;
  q->schema.effects = q->effects; q->schema.effect_count = ne;
}

/* ---- breeder's equation --------------------------------------------------- */

static void test_breeder(void) {
  quant_schema *q = (quant_schema *)malloc(sizeof *q);
  fbs_gen_config cfg = fbs_gen_config_default();
  fbs_gen_id *ids = (fbs_gen_id *)malloc(sizeof *ids * 4000u), *kids = (fbs_gen_id *)malloc(sizeof *kids * 4000u);
  fbs_gen_scored *top = (fbs_gen_scored *)malloc(sizeof *top * 4000u);
  fbs_gen_contributor *con = (fbs_gen_contributor *)malloc(sizeof *con * 800u);
  fbs_gen_union *un = (fbs_gen_union *)malloc(sizeof *un * 4000u);
  double sum_ratio = 0.0;
  uint32_t seed, i, nseeds = 5;
  if (!q || !ids || !kids || !top || !con || !un) { CHECK(0, "alloc"); goto out; }
  quant_build(q, 20, 30, 268, 128);          /* V_G = 20 * 4 * 30^2 = 72,000 = 268^2: h^2 = 0.5 */
  cfg.max_living = 8000; cfg.max_records = 8000; cfg.max_brood = 4000; cfg.max_contributors = 800;
  for (seed = 0; seed < nseeds; seed++) {
    void *mem;
    fbs_gen *g = make_ctx(&cfg, &q->schema, &mem);
    fbs_gen_index ix;
    fbs_gen_trait_stats all, sel, kid;
    fbs_gen_pool_desc pd;
    uint32_t written;
    size_t n;
    double h2, S, R, pred;
    if (!g) { CHECK(0, "breeder ctx"); continue; }
    CHECK_OK(fbs_gen_add_founders(g, 4000, 0, 100u + seed, ids));
    memset(&ix, 0, sizeof ix); ix.weight[0] = 1000;
    CHECK_OK(fbs_gen_rank(g, &ix, NULL, 0, top, 800, &n));
    for (i = 0; i < 800u; i++) { con[i].id = top[i].id; con[i].weight = 1; }
    CHECK_OK(fbs_gen_stats(g, 0, NULL, 0, &all));
    CHECK_OK(fbs_gen_stats(g, 0, ids, 4000, &sel));
    CHECK(memcmp(&sel, &all, sizeof sel) == 0, "subset of everyone equals everyone");
    {
      fbs_gen_id sid[800];
      for (i = 0; i < 800u; i++) sid[i] = top[i].id;
      CHECK_OK(fbs_gen_stats(g, 0, sid, 800, &sel));
    }
    memset(&pd, 0, sizeof pd);
    pd.children = 4000; pd.seed = 500u + seed; pd.kin_max_q16 = 65536; pd.redraws = 16;
    CHECK_OK(fbs_gen_plan_pool(g, &pd, con, 800, un, 4000, &written));
    {
      fbs_gen_brood_desc bd;
      memset(&bd, 0, sizeof bd); bd.seed = 900u + seed;
      CHECK_OK(fbs_gen_brood(g, &bd, un, written, kids, NULL));
    }
    CHECK_OK(fbs_gen_stats(g, 0, kids, written, &kid));
    h2 = (double)all.h2_q16 / 65536.0;
    S = (double)(sel.mean - all.mean);
    R = (double)(kid.mean - all.mean);
    pred = h2 * S;
    sum_ratio += R / pred;
    printf("breeder seed %u: h2 %.3f S %.1f predicted R %.1f realized R %.1f (ratio %.3f)\n",
           100u + seed, h2, S, pred, R, R / pred);
    CHECK(h2 > 0.40 && h2 < 0.60, "founder h2 %.3f near 0.5", h2);
    CHECK(fabs(R - pred) < 0.12 * pred, "seed %u realized %.1f vs predicted %.1f", seed, R, pred);
    free(mem);
  }
  printf("breeder mean realized/predicted over %u seeds: %.3f\n", nseeds, sum_ratio / nseeds);
  CHECK(fabs(sum_ratio / nseeds - 1.0) < 0.06, "mean ratio %.3f", sum_ratio / nseeds);
out:
  free(q); free(ids); free(kids); free(top); free(con); free(un);
}

/* ---- forecast vs committed children --------------------------------------- */

static void forecast_check(fbs_gen *g, fbs_gen_id a, fbs_gen_id b, uint32_t T, uint32_t N, int strict,
                           const char *label) {
  fbs_gen_forecast fc[GOB_TRAITS > 2 ? GOB_TRAITS : 2], fp[GOB_TRAITS > 2 ? GOB_TRAITS : 2];
  fbs_gen_union *un = (fbs_gen_union *)malloc(sizeof *un * N);
  fbs_gen_id *kids = (fbs_gen_id *)malloc(sizeof *kids * N);
  fbs_gen_brood_desc bd;
  uint32_t i, t, f;
  if (!un || !kids) { CHECK(0, "alloc"); free(un); free(kids); return; }
  CHECK_OK(fbs_gen_forecast_union(g, a, b, NULL, fc, T, &f));
  for (i = 0; i < N; i++) { un[i].a = a; un[i].b = b; un[i].stream = i; un[i].kinship_q30 = f; }
  CHECK_OK(fbs_gen_forecast_plan(g, un, 1, NULL, fp, T));
  CHECK(memcmp(fc, fp, sizeof fc[0] * T) == 0, "%s: forecast_plan of one union equals forecast_union", label);
  memset(&bd, 0, sizeof bd); bd.seed = 4242;
  bd.env = NULL;
  {
    fbs_gen_environment env;
    fbs_gen_environment_default(&env);
    env.mutation_scale_q8 = 0;   /* the forecast ignores mutation */
    bd.env = &env;
    CHECK_OK(fbs_gen_brood(g, &bd, un, N, kids, NULL));
  }
  for (t = 0; t < T; t++) {
    fbs_gen_trait_stats st;
    double se, derr;
    CHECK_OK(fbs_gen_stats(g, t, kids, N, &st));
    se = (double)st.sd / sqrt((double)N);
    derr = st.sd ? ((double)fc[t].sd - (double)st.sd) / (double)st.sd : 0.0;
    printf("forecast %s trait %u: forecast mean %d sd %d | children mean %d sd %d | mean %+.2f SE, sd %+.1f%%\n",
           label, t, fc[t].mean, fc[t].sd, st.mean, st.sd, se > 0 ? ((double)st.mean - fc[t].mean) / se : 0.0,
           100.0 * derr);
    if (strict) {
      CHECK(fabs((double)st.mean - (double)fc[t].mean) <= 3.0 * se + 0.5, "%s trait %u mean", label, t);
      CHECK(fabs(derr) < 0.05, "%s trait %u sd %d vs %d", label, t, fc[t].sd, st.sd);
    } else {
      CHECK(fabs((double)st.mean - (double)fc[t].mean) <= 5.0 * se + 5.0, "%s trait %u mean (loose)", label, t);
      CHECK(fabs(derr) < 0.25, "%s trait %u sd %d vs %d (loose)", label, t, fc[t].sd, st.sd);
    }
    CHECK(fc[t].p10 <= fc[t].mean && fc[t].mean <= fc[t].p90, "p10 <= mean <= p90");
  }
  free(un); free(kids);
}

static void test_forecast(void) {
  quant_schema *q = (quant_schema *)malloc(sizeof *q);
  goblin_schema_data *gs = (goblin_schema_data *)malloc(sizeof *gs);
  fbs_gen_config cfg = fbs_gen_config_default();
  fbs_gen_id ids[64];
  void *mem;
  fbs_gen *g;
  uint32_t seed;
  if (!q || !gs) { CHECK(0, "alloc"); free(q); free(gs); return; }
  quant_build(q, 30, 25, 150, 220);          /* half the loci partly dominant on trait 1 */
  cfg.max_living = 20100; cfg.max_records = 20100; cfg.max_brood = 20000; cfg.max_contributors = 64;
  for (seed = 1; seed <= 2u; seed++) {
    g = make_ctx(&cfg, &q->schema, &mem);
    if (!g) { CHECK(0, "forecast ctx"); continue; }
    CHECK_OK(fbs_gen_add_founders(g, 2, 0, seed, ids));
    forecast_check(g, ids[0], ids[1], 2, 20000, 1, "unlinked");
    free(mem);
  }
  goblin_schema_build(gs, 400);               /* common giant-kin alleles so the major loci segregate */
  g = make_ctx(&cfg, &gs->schema, &mem);
  if (g) {
    CHECK_OK(fbs_gen_add_founders(g, 64, 0, 3, ids));
    forecast_check(g, ids[0], ids[1], GOB_TRAITS, 20000, 0, "goblin");
    free(mem);
  }
  free(q); free(gs);
}

/* ---- troll line against the prototype -------------------------------------- */

static void test_troll_line(void) {
  goblin_run_config rc = goblin_run_config_a();
  goblin_run_result *res = (goblin_run_result *)malloc(sizeof *res);
  uint64_t seed;
  uint32_t s;
  if (!res) { CHECK(0, "alloc"); return; }
  for (seed = 1; seed <= 3u; seed++) {
    const goblin_row *last;
    CHECK(goblin_run_troll(seed, &rc, res) == FBS_GEN_OK, "troll run seed %u: %s", (unsigned)seed,
          fbs_gen_status_name(res->status));
    printf("troll line A, seed %u\nseason  line brood  size strength deft stride climb stamina fec     F  troll  g1   g2  wdeft\n",
           (unsigned)seed);
    for (s = 0; s < rc.seasons; s++) {
      const goblin_row *r = &res->rows[s];
      printf("%6u %5u %5u %5.2f %8.2f %4.2f %6.2f %5.2f %7.2f %4.2f %5.3f %4.0f%% %4.2f %4.2f %5.2f\n",
             r->season, r->line_n, r->brood, r->size / 1000.0, r->strength / 1000.0, r->deftness / 1000.0,
             r->stride / 1000.0, r->climb / 1000.0, r->stamina / 1000.0, r->fecundity / 1000.0,
             (double)r->f_q30 / 1073741824.0, r->troll_permille / 10.0, r->giant1_q16 / 65536.0,
             r->giant2_q16 / 65536.0, r->worker_deftness / 1000.0);
    }
    last = &res->rows[rc.seasons - 1u];
    printf("first season with line mean size >= 2.5: %u\n", res->first_season_2500);
    /* prototype (5 seeds): first season 8-10; season 12 size 2.93, F 0.076, fecundity 0.56,
     * deftness 0.62, share >= 2.5 92%, worker deftness 1.35 */
    CHECK(res->first_season_2500 >= 6u && res->first_season_2500 <= 12u, "first season %u", res->first_season_2500);
    CHECK(last->size >= 2500 && last->size <= 3300, "season-12 size %d", last->size);
    CHECK(last->f_q30 > 1073741824u / 100u && last->f_q30 < 1073741824u / 7u, "season-12 F %u", last->f_q30);
    CHECK(last->fecundity < 750, "season-12 fecundity %d", last->fecundity);
    CHECK(last->deftness < 800, "season-12 deftness %d", last->deftness);
    CHECK(last->troll_permille >= 600u, "season-12 troll share %u", last->troll_permille);
    CHECK(last->worker_deftness > 1150, "worker deftness %d", last->worker_deftness);
    CHECK(last->line_n >= 60u, "line size %u", last->line_n);
  }
  free(res);
}

int main(void) {
  test_mendel();
  test_linkage();
  test_breeder();
  test_forecast();
  test_troll_line();
  printf("%d checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
