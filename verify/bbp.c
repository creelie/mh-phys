/* Parallel branch and bound proof of F_R(g) >= T for all g in [0,inf)^(R-1), where
   F_R(g) = mu0*sum g + sum_{r=1}^{R-1} a_r sum_i w(y_{i+r}-y_i),  a_r = 2/(R-r),  y_0 = 0, y_i = g_1+...+g_i.
   Input: table of rigorous lower bounds of w on cells [c/G,(c+1)/G] (table.py, Arb at 128 bits).
   Every thread computes in round-toward-minus-infinity mode on nonnegative quantities, so every bound is a lower bound;
   T is rounded up.  A box of cells is discarded when its lower bound reaches T, or when the sum of its gap lower
   ends reaches the cutoff C with mu0*C/G >= T.  The proof fails if a single-cell box has lower bound below T.
   The work boxes after the presplit are deterministic; with a shard (index, count) only the boxes j with
   j mod count = index are searched, and every shard prints the same box count and checksum of the box list,
   so the shards together cover the whole list.
   usage: bbp R G table.bin mu0_num mu0_den T_num T_den [presplit_target [shard_index shard_count]] */
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <omp.h>
#pragma STDC FENV_ACCESS ON
#define MAXR 16
typedef struct { int lo[MAXR], hi[MAXR]; } Box;
static int R, G, n; static long L, C; static double *tab, **sp, T, mug, a[MAXR];
static double rmq(long l, long r) {
  if (l >= L || r >= L) return 0.0;                     /* beyond the table only w >= 0 is used */
  long len = r - l + 1; int k = 63 - __builtin_clzl(len); double x = sp[k][l], y = sp[k][r - (1L<<k) + 1]; return x < y ? x : y; }
/* lower bound of F_R on a box; returns 1e300 if the linear term alone settles it */
static double lower(const Box *b) {
  long slo = 0; for (int i = 0; i < n; i++) slo += b->lo[i];
  if (slo >= C) return 1e300;
  double lb = mug * (double)slo;
  long plo[MAXR+1], phi[MAXR+1]; plo[0] = phi[0] = 0;
  for (int i = 0; i < n; i++) { plo[i+1] = plo[i] + b->lo[i]; phi[i+1] = phi[i] + b->hi[i] + 1; }
  for (int r = 1; r < R && lb < T; r++) for (int i = 0; i + r <= n && lb < T; i++)
    lb += a[r] * rmq(plo[i+r] - plo[i], phi[i+r] - phi[i] - 1);
  return lb; }
static int widest(const Box *b) { int w = -1, bw = 0; for (int i = 0; i < n; i++) if (b->hi[i] - b->lo[i] > bw) { bw = b->hi[i] - b->lo[i]; w = i; } return w; }
int main(int argc, char **argv) {
  if (argc < 8) { fprintf(stderr, "usage\n"); return 2; }
  fesetround(FE_DOWNWARD);
  R = atoi(argv[1]); G = atoi(argv[2]); n = R - 1;
  FILE *f = fopen(argv[3], "rb"); if (!f) return 3; fseek(f, 0, SEEK_END); L = ftell(f)/8; fseek(f, 0, SEEK_SET);
  tab = malloc(L*8); if (fread(tab, 8, L, f) != (size_t)L) return 3; fclose(f);
  int levels = 64 - __builtin_clzl(L); sp = malloc(levels*sizeof(double*)); sp[0] = tab;
  for (int k = 1; k < levels; k++) { long w = 1L<<k; sp[k] = malloc(L*8); for (long i = 0; i + w <= L; i++) { double x = sp[k-1][i], y = sp[k-1][i + w/2]; sp[k][i] = x < y ? x : y; } }
  long mn = atol(argv[4]), md = atol(argv[5]), tn = atol(argv[6]), td = atol(argv[7]);
  long target = argc > 8 ? atol(argv[8]) : 200000;
  long sidx = argc > 10 ? atol(argv[9]) : 0, scnt = argc > 10 ? atol(argv[10]) : 1;
  if (scnt < 1 || sidx < 0 || sidx >= scnt) { fprintf(stderr, "bad shard\n"); return 2; }
  fesetround(FE_UPWARD); T = (double)tn / (double)td; fesetround(FE_DOWNWARD);
  mug = ((double)mn / (double)md) / (double)G;
  for (int r = 1; r < R; r++) a[r] = 2.0 / (double)(R - r);
  C = 0; while (mug * (double)C < T) C++;
  int ncomp = 0, cl[256], ch[256];
  for (long c = 0; c <= C; c++) { double u = mug*(double)c + a[1]*(c < L ? tab[c] : 0.0);
    if (u < T) { if (ncomp && ch[ncomp-1] == c - 1) ch[ncomp-1] = c; else { cl[ncomp] = c; ch[ncomp] = c; ncomp++; } } }
  printf("R=%d G=%d cells=%ld T=%.17g cutoff_cells=%ld components=%d:", R, G, L, T, C, ncomp);
  for (int i = 0; i < ncomp; i++) printf(" [%d,%d]", cl[i], ch[i]); printf("\n"); fflush(stdout);
  /* initial boxes: products of components, then breadth-first splitting to about `target` live boxes */
  long cap = 1L<<24; Box *q = malloc(cap*sizeof(Box)); long nq = 0;
  int idx[MAXR] = {0};
  for (;;) { Box b; for (int i = 0; i < n; i++) { b.lo[i] = cl[idx[i]]; b.hi[i] = ch[idx[i]]; } if (lower(&b) < T) q[nq++] = b;
    int i = 0; while (i < n && ++idx[i] == ncomp) { idx[i] = 0; i++; } if (i == n) break; }
  long presplit_nodes = 0;
  while (nq > 0 && nq < target) {
    long m = nq; Box *q2 = malloc(cap*sizeof(Box)); long n2 = 0;
    for (long j = 0; j < m; j++) { Box b = q[j]; int w = widest(&b);
      if (w < 0) { q2[n2++] = b; continue; }
      int mid = (b.lo[w] + b.hi[w]) / 2; Box c = b; b.hi[w] = mid; c.lo[w] = mid + 1; presplit_nodes += 2;
      if (n2 + 2 > cap) { fprintf(stderr, "queue overflow\n"); return 5; }
      if (lower(&b) < T) q2[n2++] = b; if (lower(&c) < T) q2[n2++] = c; }
    free(q); q = q2; if (n2 == m) break; nq = n2; }
  unsigned long long hsum = 1469598103934665603ULL;   /* FNV-1a over the box list */
  for (long j = 0; j < nq; j++) for (int i = 0; i < n; i++) { hsum = (hsum ^ (unsigned)q[j].lo[i]) * 1099511628211ULL; hsum = (hsum ^ (unsigned)q[j].hi[i]) * 1099511628211ULL; }
  printf("work boxes after presplit: %ld checksum=%016llx\n", nq, hsum);
  long nmine = 0; for (long j = sidx; j < nq; j += scnt) q[nmine++] = q[j];
  if (scnt > 1) printf("shard %ld/%ld: %ld work boxes\n", sidx, scnt, nmine);
  fflush(stdout);
  /* Work sharing: a global pool holds boxes still to be searched.  Every thread searches depth first from a
     box taken from the pool, and when some thread is idle it moves the older half of its stack (the larger boxes)
     back to the pool.  The tree of boxes does not depend on the order, so the node count is deterministic. */
  long total_nodes = presplit_nodes; volatile int failed = 0; volatile int nidle = 0; long shared_nodes = 0, next_report = 1000000000L;
  long gcap = nmine > 1024 ? nmine : 1024, gn = nmine; Box *pool = malloc(gcap*sizeof(Box));
  for (long j = 0; j < nmine; j++) pool[j] = q[nmine - 1 - j];      /* popped from the end: list order */
  omp_lock_t lock; omp_init_lock(&lock);
  #pragma omp parallel reduction(+:total_nodes)
  { fesetround(FE_DOWNWARD);
    int nth = omp_get_num_threads();
    long scap = 1L<<20, top = 0, since = 0; Box *st = malloc(scap*sizeof(Box));
    for (;;) {
      int got = 0, quit = 0, counted = 0;
      while (!got && !quit) {
        omp_set_lock(&lock);
        if (failed) quit = 1;
        else if (gn > 0) { st[top++] = pool[--gn]; got = 1; if (counted) nidle--; }
        else { if (!counted) { nidle++; counted = 1; } if (nidle == nth) quit = 1; }
        omp_unset_lock(&lock); }
      if (quit) break;
      while (top > 0 && !failed) {
        Box b = st[--top]; total_nodes++;
        if (++since == (1L<<20)) { since = 0; long s2;
          #pragma omp atomic capture
          { shared_nodes += 1L<<20; s2 = shared_nodes; }
          if (s2 >= next_report) {
            #pragma omp critical
            { if (s2 >= next_report) { printf("progress %ld nodes\n", s2); fflush(stdout); next_report += 1000000000L; } } }
          int idle;
          #pragma omp atomic read
          idle = nidle;
          if (idle > 0 && top >= 2) {
            omp_set_lock(&lock); long h = top/2;
            if (gn + h > gcap) { while (gn + h > gcap) gcap *= 2; pool = realloc(pool, gcap*sizeof(Box)); }
            for (long j = 0; j < h; j++) pool[gn++] = st[j];
            for (long j = h; j < top; j++) st[j - h] = st[j]; top -= h;
            omp_unset_lock(&lock); } }
        double lb = lower(&b); if (lb >= T) continue;
        int w = widest(&b);
        if (w < 0) {
          #pragma omp critical
          { if (!failed) { failed = 1; printf("FAIL box:"); for (int i = 0; i < n; i++) printf(" %d", b.lo[i]); printf(" lb=%.10g\n", lb); fflush(stdout); } }
          break; }
        int mid = (b.lo[w] + b.hi[w]) / 2; Box c = b; b.hi[w] = mid; c.lo[w] = mid + 1;
        if (top + 2 >= scap) { scap *= 2; st = realloc(st, scap*sizeof(Box)); }
        st[top++] = b; st[top++] = c; }
      if (failed) break; }
    free(st); }
  omp_destroy_lock(&lock);
  if (failed) return 1;
  if (scnt > 1) printf("VERIFIED SHARD %ld/%ld R=%d: F >= %ld/%ld with mu0=%ld/%ld on %ld of %ld work boxes (checksum=%016llx); grid=%d threads=%d nodes=%ld\n",
                       sidx, scnt, R, tn, td, mn, md, nmine, nq, hsum, G, omp_get_max_threads(), total_nodes);
  else printf("VERIFIED R=%d: F >= %ld/%ld with mu0=%ld/%ld; grid=%d threads=%d nodes=%ld\n", R, tn, td, mn, md, G, omp_get_max_threads(), total_nodes);
  return 0; }
