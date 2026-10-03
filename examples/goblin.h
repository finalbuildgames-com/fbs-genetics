/* examples/goblin.h — the goblin trait preset (design.md section 7) as header-only data.
 *
 * Host data, not part of the library. C99 and C++ compatible, no allocation. A host includes it
 * once, calls goblin_schema_build() and passes `schema` to fbs_gen_memory_for / fbs_gen_create;
 * the module copies the schema, so the goblin_schema_data may be discarded afterwards.
 * Shared by the tests and examples/gathering.c so they run the same preset.
 *
 * Numbers match Studio/Research/Genetics-2026-09-26/prototype/goblin_breeding_sim.py: 42 trait
 * loci on 6 chromosomes of 100 cM (map order and Haldane recombination fractions computed from the
 * prototype's layout), plus the six cosmetic loci of section 7.2 placed at each chromosome's end.
 * env_sd comes from the founder genetic variance and the target h^2 (AlphaSimR's rule), computed
 * analytically: quantitative locus var = 4 w^2 (uniform -2..2 alleles), major locus var =
 * 2 p (1 - p) w^2 with p = 0.04.
 */
#ifndef FBS_GENETICS_EXAMPLES_GOBLIN_H
#define FBS_GENETICS_EXAMPLES_GOBLIN_H

#include "fbs/genetics.h"

#include <string.h>

enum { GOB_SIZE, GOB_STRENGTH, GOB_DEFTNESS, GOB_STRIDE, GOB_STAMINA, GOB_THRIFT, GOB_CLIMB,
       GOB_HEAT, GOB_DAMP, GOB_FECUNDITY, GOB_TRAITS };

static const fbs_gen_trait_desc goblin_traits[GOB_TRAITS] = {
  /* key            len  base   lo    hi   env_sd depression */
  {"size",           4, 1000,  600, 4000, 232,   0},
  {"strength",       8, 1000,  400, 5000, 294,   0},
  {"deftness",       8, 1000,  400, 1800, 165, 200},
  {"stride",         6, 1000,  500, 1600, 147,   0},
  {"stamina",        7, 1000,  400, 2000, 230, 400},
  {"thrift",         6, 1000,  500, 1600, 180,   0},
  {"climb",          5, 1000,  200, 1800, 258,   0},
  {"heat",           4,  500,    0, 1000, 294,   0},
  {"damp",           4,  500,    0, 1000, 294,   0},
  {"fecundity",      9, 1000,  100, 2000, 190, 500},
};

/* Effect groups: every locus of a group carries the same effects. */
enum { G_SZ, G_GIANT, G_ST, G_DF, G_RN, G_EN, G_TH, G_CL, G_HEATP, G_DAMPP, G_HT, G_DM, G_FC,
       G_EAR, G_MOTTLE, G_EYEGLOW, G_TUSK, G_GROUPS };

typedef struct goblin_group_effect { uint8_t group, trait; int32_t weight; } goblin_group_effect;
static const goblin_group_effect goblin_group_effects[] = {
  {G_SZ, GOB_SIZE, 30}, {G_SZ, GOB_DEFTNESS, -6}, {G_SZ, GOB_CLIMB, -10}, {G_SZ, GOB_FECUNDITY, -4},
  {G_SZ, GOB_STRIDE, -3},
  {G_GIANT, GOB_SIZE, 450}, {G_GIANT, GOB_STRENGTH, 350}, {G_GIANT, GOB_DEFTNESS, -80},
  {G_GIANT, GOB_STRIDE, -50}, {G_GIANT, GOB_CLIMB, -150}, {G_GIANT, GOB_FECUNDITY, -120},
  {G_GIANT, GOB_STAMINA, 50},
  {G_ST, GOB_STRENGTH, 40}, {G_ST, GOB_THRIFT, -8}, {G_ST, GOB_DEFTNESS, -4},
  {G_DF, GOB_DEFTNESS, 25}, {G_DF, GOB_STRENGTH, -5},
  {G_RN, GOB_STRIDE, 25}, {G_RN, GOB_STAMINA, -10},
  {G_EN, GOB_STAMINA, 40}, {G_EN, GOB_THRIFT, -10},
  {G_TH, GOB_THRIFT, 30}, {G_TH, GOB_STRIDE, -8},
  {G_CL, GOB_CLIMB, 50}, {G_CL, GOB_SIZE, -10},
  {G_HEATP, GOB_HEAT, 50}, {G_HEATP, GOB_DAMP, -50},
  {G_DAMPP, GOB_DAMP, 50}, {G_DAMPP, GOB_HEAT, -50},
  {G_HT, GOB_HEAT, 40}, {G_DM, GOB_DAMP, 40},
  {G_FC, GOB_FECUNDITY, 40}, {G_FC, GOB_STAMINA, -10},
};

typedef struct goblin_map_entry { const char *key; uint8_t group; uint16_t rec; } goblin_map_entry;
#define GOB_LOCI 48
static const goblin_map_entry goblin_map[GOB_LOCI] = {
  /* chromosome 1 */
  {"SZ1", G_SZ, 32768}, {"GIANT1", G_GIANT, 6529}, {"ST6", G_ST, 6529}, {"DF5", G_DF, 6529},
  {"RN4", G_RN, 6529}, {"EN3", G_EN, 6529}, {"TH2", G_TH, 6529}, {"CL1", G_CL, 6529},
  {"HT1", G_HT, 6529}, {"EAR", G_EAR, 3446},
  /* chromosome 2 */
  {"SZ2", G_SZ, 32768}, {"ST1", G_ST, 8144}, {"EN4", G_EN, 8144}, {"TH3", G_TH, 8144},
  {"CL2", G_CL, 8144}, {"HD1", G_HEATP, 8144}, {"FC1", G_FC, 8144}, {"SNOUT", G_EAR, 4362},
  /* chromosome 3 */
  {"SZ3", G_SZ, 32768}, {"ST2", G_ST, 8144}, {"DF1", G_DF, 8144}, {"TH4", G_TH, 8144},
  {"CL3", G_CL, 8144}, {"HD3", G_DAMPP, 8144}, {"FC2", G_FC, 8144}, {"HEAD", G_EAR, 4362},
  /* chromosome 4 */
  {"SZ4", G_SZ, 32768}, {"GIANT2", G_GIANT, 9289}, {"ST3", G_ST, 9289}, {"DF2", G_DF, 9289},
  {"RN1", G_RN, 9289}, {"DM1", G_DM, 9289}, {"MOTTLE", G_MOTTLE, 5030},
  /* chromosome 5 */
  {"SZ5", G_SZ, 32768}, {"ST4", G_ST, 9289}, {"DF3", G_DF, 9289}, {"RN2", G_RN, 9289},
  {"EN1", G_EN, 9289}, {"HD2", G_HEATP, 9289}, {"EYEGLOW", G_EYEGLOW, 5030},
  /* chromosome 6 */
  {"SZ6", G_SZ, 32768}, {"ST5", G_ST, 8144}, {"DF4", G_DF, 8144}, {"RN3", G_RN, 8144},
  {"EN2", G_EN, 8144}, {"TH1", G_TH, 8144}, {"HD4", G_DAMPP, 8144},
  {"TUSK", G_TUSK, 4362},
};

typedef struct goblin_schema_data {
  fbs_gen_locus_desc loci[GOB_LOCI];
  fbs_gen_effect effects[GOB_LOCI * 8];
  fbs_gen_schema schema;
} goblin_schema_data;

/* giant_permille: founder frequency of the giant-kin allele (40 = 4%). */
static void goblin_schema_build(goblin_schema_data *d, uint16_t giant_permille) {
  uint32_t i, j, ne = 0;
  memset(d, 0, sizeof *d);
  for (i = 0; i < GOB_LOCI; i++) {
    const goblin_map_entry *m = &goblin_map[i];
    fbs_gen_locus_desc *l = &d->loci[i];
    l->key = m->key; l->key_len = strlen(m->key);
    l->recombination_q16 = m->rec;
    l->h = 128;
    switch (m->group) {
      case G_GIANT:
        l->amin = 0; l->amax = 1;
        l->founder_count = 2;
        l->founder[0].allele = 0; l->founder[0].weight = (uint16_t)(1000u - giant_permille);
        l->founder[1].allele = 1; l->founder[1].weight = giant_permille;
        break;
      case G_MOTTLE:
        l->amin = 0; l->amax = 2; l->founder_count = 3;
        l->founder[0].allele = 0; l->founder[0].weight = 70;
        l->founder[1].allele = 1; l->founder[1].weight = 20;
        l->founder[2].allele = 2; l->founder[2].weight = 10;
        break;
      case G_EYEGLOW: /* recessive glow: the lower allele dominates */
        l->amin = 0; l->amax = 1; l->h = 0; l->founder_count = 2;
        l->founder[0].allele = 0; l->founder[0].weight = 90;
        l->founder[1].allele = 1; l->founder[1].weight = 10;
        break;
      case G_TUSK:
        l->amin = 0; l->amax = 1; l->h = 256; l->founder_count = 2;
        l->founder[0].allele = 0; l->founder[0].weight = 80;
        l->founder[1].allele = 1; l->founder[1].weight = 20;
        break;
      default: /* quantitative (and the ear/snout/head shape loci) */
        l->amin = (int8_t)(m->group == G_EAR ? -4 : -8);
        l->amax = (int8_t)(m->group == G_EAR ? 4 : 8);
        l->mutation_q16 = 66;   /* 0.1% per allele per gamete */
        l->mutation_step = 1;
        l->founder_count = 5;
        for (j = 0; j < 5; j++) { l->founder[j].allele = (int8_t)((int)j - 2); l->founder[j].weight = 1; }
        break;
    }
    for (j = 0; j < sizeof goblin_group_effects / sizeof goblin_group_effects[0]; j++) {
      const goblin_group_effect *g = &goblin_group_effects[j];
      if (g->group != m->group) continue;
      d->effects[ne].locus = (uint16_t)i;
      d->effects[ne].trait = g->trait;
      d->effects[ne].match = m->group == G_GIANT ? (int8_t)1 : (int8_t)FBS_GEN_ANY;
      d->effects[ne].weight = g->weight;
      ne++;
    }
  }
  d->schema.traits = goblin_traits; d->schema.trait_count = GOB_TRAITS;
  d->schema.loci = d->loci; d->schema.locus_count = GOB_LOCI;
  d->schema.effects = d->effects; d->schema.effect_count = ne;
}

#endif /* FBS_GENETICS_EXAMPLES_GOBLIN_H */
