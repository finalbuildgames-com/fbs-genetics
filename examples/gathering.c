/* examples/gathering.c — human-readable dump of a 12-season troll-line run on the goblin preset.
 *
 * Prints the line means per season in the shape of design.md section 8 and a lineage tablet (sire
 * and dam to depth 3) for the biggest living troll, so a person can read what the assertions in
 * tests/ did not predict. Host code: it allocates.
 *
 * Usage: fbs-genetics-gathering [seed]     (default seed 1)
 */
#include "fbs/genetics.h"
#include "goblin_troll.h"

#include <stdio.h>
#include <stdlib.h>

static const char *kin_mark(uint32_t f_q30) {
  if (f_q30 >= FBS_GEN_Q30 / 4u) return "<<<";
  if (f_q30 >= FBS_GEN_Q30 / 8u) return "<<";
  if (f_q30 >= FBS_GEN_Q30 / 16u) return "<";
  return "";
}

static void tablet_line(fbs_gen *g, const char *role, fbs_gen_id id, int indent) {
  fbs_gen_info info;
  int32_t z[GOB_TRAITS];
  if (id == 0u) { printf("%*s%-10s unknown (founder stock)\n", indent, "", role); return; }
  if (fbs_gen_get(g, id, &info) != FBS_GEN_OK) {
    printf("%*s%-10s #%u (record forgotten)\n", indent, "", role, id);
    return;
  }
  printf("%*s%-10s #%-5u %-6s line %u  born %2u  ", indent, "", role, id,
         (info.flags & FBS_GEN_ALIVE) ? "living" : "dead", info.line, info.born);
  if (fbs_gen_birth_traits(g, id, z, GOB_TRAITS) == FBS_GEN_OK)
    printf("size %.2f strength %.2f deft %.2f fec %.2f  ", z[GOB_SIZE] / 1000.0, z[GOB_STRENGTH] / 1000.0,
           z[GOB_DEFTNESS] / 1000.0, z[GOB_FECUNDITY] / 1000.0);
  printf("F %.3f %-3s children %u honoured %u\n", (double)info.f_q30 / FBS_GEN_Q30, kin_mark(info.f_q30),
         info.children, info.honours);
}

int main(int argc, char **argv) {
  goblin_run_config rc = goblin_run_config_a();
  goblin_run_result *res = (goblin_run_result *)malloc(sizeof *res);
  uint64_t seed = argc > 1 ? (uint64_t)strtoull(argv[1], NULL, 10) : 1u;
  fbs_gen *g = NULL;
  void *mem = NULL;
  uint32_t s;
  int st;
  if (!res) return 1;
  st = goblin_run_troll_keep(seed, &rc, res, &g, &mem);
  if (st != FBS_GEN_OK) {
    fprintf(stderr, "troll run failed: %s\n", fbs_gen_status_name(st));
    free(res);
    return 1;
  }
  printf("Troll line, config A (honour 20 of the line, kin redraw above 0.125 up to 4 tries, brood 60 x\n"
         "mean fecundity; workers: 200 honoured of 1,000 by deftness + stamina, 400 children), seed %u\n\n",
         (unsigned)seed);
  printf("season  line  brood  size  strength  deftness  stride  climb  stamina  fecundity      F  size>=2.5"
         "  giant1  giant2  worker deft\n");
  for (s = 0; s < rc.seasons; s++) {
    const goblin_row *r = &res->rows[s];
    printf("%6u %5u %6u %5.2f %9.2f %9.2f %7.2f %6.2f %8.2f %10.2f %6.3f %9.1f%% %7.2f %7.2f %12.2f\n",
           r->season, r->line_n, r->brood, r->size / 1000.0, r->strength / 1000.0, r->deftness / 1000.0,
           r->stride / 1000.0, r->climb / 1000.0, r->stamina / 1000.0, r->fecundity / 1000.0,
           (double)r->f_q30 / FBS_GEN_Q30, r->troll_permille / 10.0, r->giant1_q16 / 65536.0,
           r->giant2_q16 / 65536.0, r->worker_deftness / 1000.0);
  }
  if (res->first_season_2500) printf("\nline mean size first reaches 2.5 in season %u\n", res->first_season_2500);
  else printf("\nline mean size never reaches 2.5\n");
  printf("final snapshot: %lu bytes, FNV-1a 64 %016llx\n", (unsigned long)res->snapshot_bytes,
         (unsigned long long)res->snapshot_digest);

  { /* lineage tablet of the biggest living troll */
    const fbs_gen_id *ids; const int32_t *size;
    uint32_t living, i, best = 0;
    int32_t best_size = 0;
    fbs_gen_id anc[14];
    size_t n = 0;
    static const char *roles[14] = { "sire", "dam", "sire's sire", "sire's dam", "dam's sire", "dam's dam",
      "sss", "ssd", "sds", "sdd", "dss", "dsd", "dds", "ddd" };
    fbs_gen_column(g, GOB_SIZE, &ids, &size, &living);
    for (i = 0; i < living; i++) {
      fbs_gen_info info;
      fbs_gen_get(g, ids[i], &info);
      if (info.line == GOB_TROLL_LINE && (size[i] > best_size || (size[i] == best_size && ids[i] < best))) {
        best = ids[i]; best_size = size[i];
      }
    }
    if (best) {
      printf("\nLineage tablet (depth 3), biggest living troll; birth traits, kin marks < 1/16  << 1/8  <<< 1/4\n");
      tablet_line(g, "troll", best, 0);
      fbs_gen_ancestors(g, best, 3, anc, 14, &n);
      tablet_line(g, roles[0], anc[0], 2);
      tablet_line(g, roles[2], anc[2], 4);
      tablet_line(g, "sire", anc[6], 6); tablet_line(g, "dam", anc[7], 6);
      tablet_line(g, roles[3], anc[3], 4);
      tablet_line(g, "sire", anc[8], 6); tablet_line(g, "dam", anc[9], 6);
      tablet_line(g, roles[1], anc[1], 2);
      tablet_line(g, roles[4], anc[4], 4);
      tablet_line(g, "sire", anc[10], 6); tablet_line(g, "dam", anc[11], 6);
      tablet_line(g, roles[5], anc[5], 4);
      tablet_line(g, "sire", anc[12], 6); tablet_line(g, "dam", anc[13], 6);
      {
        uint32_t k = 0;
        if (anc[0] && anc[1] && fbs_gen_kinship(g, anc[0], anc[1], &k) == FBS_GEN_OK)
          printf("kinship of sire and dam %.4f (the troll's F)\n", (double)k / FBS_GEN_Q30);
      }
    }
  }
  free(mem);
  free(res);
  return 0;
}
