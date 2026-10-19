/* Branch and bound proof of F_R(g) >= T for all g in [0,inf)^(R-1), where
   F_R(g) = mu0*sum g + sum_{r=1}^{R-1} a_r sum_i w(y_{i+r}-y_i),  a_r = 2/(R-r).
   Input: table of rigorous lower bounds of w on cells [c/G,(c+1)/G].
   All arithmetic is in round-toward-minus-infinity mode on nonnegative quantities, so every computed bound is a lower bound.
   usage: bb R G table.bin mu0_num mu0_den T_num T_den */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#pragma STDC FENV_ACCESS ON
#define MAXR 16
static int R, G; static long L; static double *tab; static double **sp; static int levels;
static double rmq(long l, long r) { /* min of tab[l..r], 0 beyond table */
  if (l >= L) return 0.0; if (r >= L) return 0.0; /* w >= 0, and beyond the table we only know w >= 0 */
  long len = r - l + 1; int k = 63 - __builtin_clzl(len); double a = sp[k][l], b = sp[k][r - (1L<<k) + 1]; return a < b ? a : b; }
typedef struct { int lo[MAXR], hi[MAXR]; } Box;
int main(int argc, char **argv) {
  if (argc < 8) { fprintf(stderr, "usage\n"); return 2; }
  R = atoi(argv[1]); G = atoi(argv[2]); int n = R - 1;
  FILE *f = fopen(argv[3], "rb"); fseek(f, 0, SEEK_END); L = ftell(f)/8; fseek(f, 0, SEEK_SET);
  tab = malloc(L*8); if (fread(tab, 8, L, f) != (size_t)L) return 3; fclose(f);
  levels = 64 - __builtin_clzl(L); sp = malloc(levels*sizeof(double*)); sp[0] = tab;
  for (int k = 1; k < levels; k++) { long w = 1L<<k; sp[k] = malloc(L*8); for (long i = 0; i + w <= L; i++) { double a = sp[k-1][i], b = sp[k-1][i + w/2]; sp[k][i] = a < b ? a : b; } }
  long mn = atol(argv[4]), md = atol(argv[5]), tn = atol(argv[6]), td = atol(argv[7]);
  fesetround(FE_UPWARD); double T = (double)tn / (double)td;           /* T rounded up */
  fesetround(FE_DOWNWARD);
  double mug = ((double)mn / (double)md) / (double)G;                    /* mu0/G rounded down */
  double a[MAXR]; for (int r = 1; r < R; r++) a[r] = 2.0 / (double)(R - r); /* rounded down */
  /* cutoff: if sum of gap lower ends >= C cells then mu0*sum >= T */
  long C = 0; while (mug * (double)C < T) C++;
  /* per-coordinate surviving cells: mu0*g + a_1*w(g) < T */
  int *surv = calloc(C + 1, sizeof(int)); int ncomp = 0; int cl[64], ch[64];
  for (long c = 0; c <= C; c++) { double u = mug*(double)c + a[1]*(c < L ? tab[c] : 0.0); surv[c] = u < T; }
  for (long c = 0; c <= C; c++) if (surv[c]) { if (ncomp && cl[ncomp-1] <= c && ch[ncomp-1] == c - 1) ch[ncomp-1] = c; else { cl[ncomp] = c; ch[ncomp] = c; ncomp++; } }
  printf("R=%d G=%d cells=%ld T=%.17g cutoff_cells=%ld components=%d:", R, G, L, T, C, ncomp);
  for (int i = 0; i < ncomp; i++) printf(" [%d,%d]", cl[i], ch[i]); printf("\n"); fflush(stdout);
  long cap = 1L<<26; Box *st = malloc(cap*sizeof(Box)); long top = 0;
  int idx[MAXR] = {0};
  for (;;) { Box b; for (int i = 0; i < n; i++) { b.lo[i] = cl[idx[i]]; b.hi[i] = ch[idx[i]]; } st[top++] = b;
    int i = 0; while (i < n && ++idx[i] == ncomp) { idx[i] = 0; i++; } if (i == n) break; }
  long init = top, nodes = 0, maxtop = 0; double worst = 1e9;
  while (top > 0) {
    Box b = st[--top]; nodes++;
    long slo = 0; for (int i = 0; i < n; i++) slo += b.lo[i];
    if (slo >= C) continue;
    double lb = mug * (double)slo;
    long plo[MAXR+1], phi[MAXR+1]; plo[0] = phi[0] = 0;
    for (int i = 0; i < n; i++) { plo[i+1] = plo[i] + b.lo[i]; phi[i+1] = phi[i] + b.hi[i] + 1; }
    for (int r = 1; r < R && lb < T; r++) for (int i = 0; i + r <= n && lb < T; i++) {
      long l = plo[i+r] - plo[i], h = phi[i+r] - phi[i] - 1; lb += a[r] * rmq(l, h); }
    if (lb >= T) continue;
    int w = -1, bw = 0; for (int i = 0; i < n; i++) if (b.hi[i] - b.lo[i] > bw) { bw = b.hi[i] - b.lo[i]; w = i; }
    if (w < 0) { if (lb < worst) worst = lb; printf("FAIL box:"); for (int i = 0; i < n; i++) printf(" %d", b.lo[i]); printf(" lb=%.10g\n", lb); return 1; }
    int mid = (b.lo[w] + b.hi[w]) / 2; Box c = b; b.hi[w] = mid; c.lo[w] = mid + 1;
    if (top + 2 >= cap) { fprintf(stderr, "stack overflow\n"); return 4; }
    st[top++] = b; st[top++] = c; if (top > maxtop) maxtop = top;
    if ((nodes & ((1L<<27)-1)) == 0) { printf("progress nodes=%ld stack=%ld\n", nodes, top); fflush(stdout); }
  }
  printf("VERIFIED R=%d: F >= %ld/%ld with mu0=%ld/%ld; initial_boxes=%ld nodes=%ld max_stack=%ld\n", R, tn, td, mn, md, init, nodes, maxtop);
  return 0;
}
