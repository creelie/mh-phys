/* Second-order version of bbp_ww: parallel branch and bound proof of F(g) >= T for all g in [0,inf)^(R-1), where
   F(g) = sum_i b_i g_i + sum_{r=1}^{R-1} sum_i a_{r,i} w(y_{i+r}-y_i),  y_0 = 0, y_i = g_1+...+g_i,
   with weights read from a file of lines "b i num den" (i = 1..R-1) and "a r i num den" (0 <= i <= R-1-r).
   Boxes, cutoffs, presplit, shards and checksum are those of bbp_ww, and the presplit uses the cell bound of bbp_ww.
   In the search a box is discarded when the larger of two lower bounds reaches T:
   (1) the cell bound of bbp_ww, from the table of lower bounds of w on cells (table.py);
   (2) a second-order bound at the grid point c of the box nearest its centre, from the tables of dtable.py:
       for each pair, d = d_c + D with D in [D-, D+] (in cells), w(d) >= w(d_c) + w'(d_c) D + (1/2) min w'' D^2,
       min w'' taken over the cells of the pair; the linear terms are collected gap by gap, with b_i, and each
       gap's coefficient interval is minimised over the corners of the box; the quadratic term is dropped when
       min w'' >= 0 and otherwise bounded by max(D-^2, D+^2).  It is skipped if a pair reaches past the table.
   All arithmetic is in round-toward-minus-infinity mode; upper bounds are formed as -((-x) op y).  T is rounded up.
   usage: bbp_cf R G table.bin dtable.bin weights.txt T_num T_den [presplit_target [shard_index shard_count]] */
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <omp.h>
#pragma STDC FENV_ACCESS ON
#define MAXR 16
typedef struct { int lo[MAXR], hi[MAXR]; } Box;
static int R, G, n; static long L, C[MAXR]; static double *tab, **sp, T, bg[MAXR], a[MAXR][MAXR];
static double *wpt, *gpl, *gph, **hmn, bgu[MAXR], au[MAXR][MAXR];   /* w and w'/G at grid points, w''/G^2 on cells, upper weights */
static inline double upmul(double x, double y) { return -((-x) * y); }
static inline double upadd(double x, double y) { return -((-x) - y); }
static double rminh(long l, long r) { long len = r - l + 1; int k = 63 - __builtin_clzl(len); double x = hmn[k][l], y = hmn[k][r - (1L<<k) + 1]; return x < y ? x : y; }
static double rmq(long l, long r) {
  if (l >= L || r >= L) return 0.0;                     /* beyond the table only w >= 0 is used */
  long len = r - l + 1; int k = 63 - __builtin_clzl(len); double x = sp[k][l], y = sp[k][r - (1L<<k) + 1]; return x < y ? x : y; }
/* the cell bound (1); returns 1e300 if the linear term alone settles it */
static double lower(const Box *b) {
  double lb = 0.0; for (int i = 0; i < n; i++) lb += bg[i] * (double)b->lo[i];
  if (lb >= T) return 1e300;
  long plo[MAXR+1], phi[MAXR+1]; plo[0] = phi[0] = 0;
  for (int i = 0; i < n; i++) { plo[i+1] = plo[i] + b->lo[i]; phi[i+1] = phi[i] + b->hi[i] + 1; }
  for (int r = 1; r < R && lb < T; r++) for (int i = 0; i + r <= n && lb < T; i++)
    lb += a[r][i] * rmq(plo[i+r] - plo[i], phi[i+r] - phi[i] - 1);
  return lb; }
/* the second-order bound (2) */
static double lower_q(const Box *b) {
  long plo[MAXR+1], phi[MAXR+1], pc[MAXR+1], el[MAXR+1], eh[MAXR+1]; int c[MAXR]; double Dl[MAXR], Dh[MAXR];
  plo[0] = phi[0] = pc[0] = el[0] = eh[0] = 0;
  for (int i = 0; i < n; i++) { c[i] = (b->lo[i] + b->hi[i] + 1) / 2; plo[i+1] = plo[i] + b->lo[i]; phi[i+1] = phi[i] + b->hi[i] + 1; pc[i+1] = pc[i] + c[i];
    el[i+1] = el[i] + (b->lo[i] - c[i]); eh[i+1] = eh[i] + (b->hi[i] + 1 - c[i]); }
  if (phi[n] > L) return -1e300;
  double f = 0.0;
  for (int i = 0; i < n; i++) { f += bg[i] * (double)c[i]; Dl[i] = bg[i]; Dh[i] = bgu[i]; }
  for (int r = 1; r < R; r++) for (int i = 0; i + r <= n; i++) {
    if (au[r][i] == 0.0) continue;
    long dc = pc[i+r] - pc[i];
    f += a[r][i] * wpt[dc];
    double m = gpl[dc], M = gph[dc];
    double lo = m >= 0 ? a[r][i] * m : au[r][i] * m, hi = M >= 0 ? upmul(au[r][i], M) : upmul(a[r][i], M);
    for (int j = i; j < i + r; j++) { Dl[j] += lo; Dh[j] = upadd(Dh[j], hi); }
    double h2 = rminh(plo[i+r] - plo[i], phi[i+r] - phi[i] - 1);
    if (h2 < 0) { double x = (double)(el[i+r] - el[i]), y = (double)(eh[i+r] - eh[i]); double q = x*x > y*y ? x*x : y*y; f += (0.5 * au[r][i] * h2) * q; } }
  for (int i = 0; i < n; i++) {
    double xl = (double)(b->lo[i] - c[i]), xh = (double)(b->hi[i] + 1 - c[i]);
    double t = Dl[i] * xl, u;
    u = Dl[i] * xh; if (u < t) t = u;
    u = Dh[i] * xl; if (u < t) t = u;
    u = Dh[i] * xh; if (u < t) t = u;
    f += t; }
  return f; }
static int widest(const Box *b) { int w = -1, bw = 0; for (int i = 0; i < n; i++) if (b->hi[i] - b->lo[i] > bw) { bw = b->hi[i] - b->lo[i]; w = i; } return w; }
int main(int argc, char **argv) {
  if (argc < 8) { fprintf(stderr, "usage\n"); return 2; }
  fesetround(FE_DOWNWARD);
  R = atoi(argv[1]); G = atoi(argv[2]); n = R - 1;
  FILE *f = fopen(argv[3], "rb"); if (!f) return 3; fseek(f, 0, SEEK_END); L = ftell(f)/8; fseek(f, 0, SEEK_SET);
  tab = malloc(L*8); if (fread(tab, 8, L, f) != (size_t)L) return 3; fclose(f);
  int levels = 64 - __builtin_clzl(L); sp = malloc(levels*sizeof(double*)); sp[0] = tab;
  for (int k = 1; k < levels; k++) { long w = 1L<<k; sp[k] = malloc(L*8); for (long i = 0; i + w <= L; i++) { double x = sp[k-1][i], y = sp[k-1][i + w/2]; sp[k][i] = x < y ? x : y; } }
  { FILE *g = fopen(argv[4], "rb"); if (!g) return 3; fseek(g, 0, SEEK_END); if (ftell(g) != (4*L + 3)*8) { fprintf(stderr, "derivative table size\n"); return 3; } fseek(g, 0, SEEK_SET);
    double *d = malloc((4*L + 3)*8); if (fread(d, 8, 4*L + 3, g) != (size_t)(4*L + 3)) return 3; fclose(g);
    wpt = d; gpl = d + (L + 1); gph = d + 2*(L + 1); hmn = malloc(levels*sizeof(double*)); hmn[0] = d + 3*(L + 1);
    for (int k = 1; k < levels; k++) { long w = 1L<<k; hmn[k] = malloc(L*8);
      for (long i = 0; i + w <= L; i++) { double x = hmn[k-1][i], y = hmn[k-1][i + w/2]; hmn[k][i] = x < y ? x : y; } } }
  long tn = atol(argv[6]), td = atol(argv[7]);
  long target = argc > 8 ? atol(argv[8]) : 200000;
  long sidx = argc > 10 ? atol(argv[9]) : 0, scnt = argc > 10 ? atol(argv[10]) : 1;
  if (scnt < 1 || sidx < 0 || sidx >= scnt) { fprintf(stderr, "bad shard\n"); return 2; }
  fesetround(FE_UPWARD); T = (double)tn / (double)td; fesetround(FE_DOWNWARD);
  for (int i = 0; i < n; i++) bg[i] = -1.0;
  for (int r = 1; r < R; r++) for (int i = 0; i + r <= n; i++) a[r][i] = -1.0;
  { FILE *wf = fopen(argv[5], "r"); if (!wf) return 3; char kind[4]; int r, i; long num, den;
    while (fscanf(wf, "%3s", kind) == 1) {
      if (kind[0] == 'b' && fscanf(wf, "%d %ld %ld", &i, &num, &den) == 3 && i >= 1 && i <= n && num > 0 && den > 0)
        { bg[i-1] = ((double)num / (double)den) / (double)G; bgu[i-1] = -((((double)(-num)) / (double)den) / (double)G); }
      else if (kind[0] == 'a' && fscanf(wf, "%d %d %ld %ld", &r, &i, &num, &den) == 4 && r >= 1 && r < R && i >= 0 && i + r <= n && num >= 0 && den > 0)
        { a[r][i] = (double)num / (double)den; au[r][i] = -(((double)(-num)) / (double)den); }
      else { fprintf(stderr, "bad weights file\n"); return 2; } }
    fclose(wf); }
  for (int i = 0; i < n; i++) if (bg[i] < 0) { fprintf(stderr, "missing b %d\n", i+1); return 2; }
  for (int r = 1; r < R; r++) for (int i = 0; i + r <= n; i++) if (a[r][i] < 0) { fprintf(stderr, "missing a %d %d\n", r, i); return 2; }
  /* gap i alone: F >= b_i g_i + a_{1,i-1} w(g_i), so only cells where this is below T can occur */
  int ncomp[MAXR], cl[MAXR][256], ch[MAXR][256];
  printf("R=%d G=%d cells=%ld T=%.17g\n", R, G, L, T);
  for (int i = 0; i < n; i++) {
    C[i] = 0; while (bg[i] * (double)C[i] < T) C[i]++;
    ncomp[i] = 0;
    for (long c = 0; c <= C[i]; c++) { double u = bg[i]*(double)c + a[1][i]*(c < L ? tab[c] : 0.0);
      if (u < T) { int k = ncomp[i]; if (k && ch[i][k-1] == c - 1) ch[i][k-1] = c; else { if (k == 256) return 6; cl[i][k] = c; ch[i][k] = c; ncomp[i]++; } } }
    printf("  gap %d: cutoff_cells=%ld components=%d:", i+1, C[i], ncomp[i]);
    for (int k = 0; k < ncomp[i]; k++) printf(" [%d,%d]", cl[i][k], ch[i][k]); printf("\n"); }
  fflush(stdout);
  /* initial boxes: products of components, then breadth-first splitting to about `target` live boxes */
  long cap = 1L<<24; Box *q = malloc(cap*sizeof(Box)); long nq = 0;
  int idx[MAXR] = {0};
  for (;;) { Box b; for (int i = 0; i < n; i++) { b.lo[i] = cl[i][idx[i]]; b.hi[i] = ch[i][idx[i]]; } if (lower(&b) < T) q[nq++] = b;
    int i = 0; while (i < n && ++idx[i] == ncomp[i]) { idx[i] = 0; i++; } if (i == n) break; }
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
        { double l2 = lower_q(&b); if (l2 > lb) lb = l2; } if (lb >= T) continue;
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
  if (scnt > 1) printf("VERIFIED SHARD %ld/%ld R=%d: F >= %ld/%ld with weights %s on %ld of %ld work boxes (checksum=%016llx); grid=%d threads=%d nodes=%ld\n",
                       sidx, scnt, R, tn, td, argv[5], nmine, nq, hsum, G, omp_get_max_threads(), total_nodes);
  else printf("VERIFIED R=%d: F >= %ld/%ld with weights %s; grid=%d threads=%d nodes=%ld\n", R, tn, td, argv[5], G, omp_get_max_threads(), total_nodes);
  return 0; }
