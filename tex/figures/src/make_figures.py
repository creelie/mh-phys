"""Figures of the paper. usage: python3 -I make_figures.py outdir
Every figure is written at 200 dpi, at most 1600 px wide, flattened onto white and reduced to a 256-colour palette."""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm",
                     "axes.spines.top": False, "axes.spines.right": False})
s2 = np.sqrt(2.0); K0 = s2*np.sin(1/s2); H = 1.5 - (np.cos(1/s2)/np.sin(1/s2))/s2
def k(x): return 0.5*(np.sinc((np.pi*x - 1/s2)/np.pi) + np.sinc((np.pi*x + 1/s2)/np.pi))/K0
def w(x): return k(x)**2
BLUE, RED, GREEN, GREY, ORANGE = "#1f4e79", "#b03a2e", "#2e7d32", "#7f7f7f", "#d4801f"

def save(fig, name):
    p = os.path.join(out, name)
    fig.savefig(p, dpi=200, facecolor="white"); plt.close(fig)
    im = Image.open(p).convert("RGB")
    if im.width > 1600: im = im.resize((1600, round(im.height*1600/im.width)), Image.LANCZOS)
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(p, optimize=True)
    print(name, im.size, os.path.getsize(p))

# Figure 1: the kernel, and the energy of three consecutive zeros
fig = plt.figure(figsize=(8, 3.4))
ax = fig.add_subplot(1, 2, 1)
x = np.linspace(0, 4, 2000)
ax.plot(x, k(x), color=BLUE, lw=1.6, label=r"$k(x)$")
ax.plot(x, w(x), color=RED, lw=1.6, label=r"$w(x)=k(x)^2$")
ax.axhline(0, color=GREY, lw=0.6)
ax.axvspan(0, 0.55, color=ORANGE, alpha=0.15)
ax.text(0.04, 0.08, "gaps\n$<0.55$", fontsize=8, color=ORANGE)
for z in (1.0572783, 2.0301, 3.0202):
    ax.plot([z], [0], "o", ms=4, color="black")
ax.set_xlabel("distance $x$ in units of the mean spacing"); ax.set_xlim(0, 4); ax.set_ylim(-0.25, 1.05)
ax.legend(frameon=False, fontsize=9); ax.set_title("(a) the overlap kernel", fontsize=10)
ax = fig.add_subplot(1, 2, 2, projection="3d")
g1, g2 = np.meshgrid(np.linspace(0.6, 3.2, 160), np.linspace(0.6, 3.2, 160))
E = np.log10(2*(w(g1) + w(g2) + w(g1 + g2)) + 1e-4)
ax.plot_surface(g1, g2, E, cmap=cm.viridis, rstride=2, cstride=2, linewidth=0, antialiased=True, alpha=0.95)
ax.set_xlabel("$g_1$", labelpad=-2); ax.set_ylabel("$g_2$", labelpad=-2); ax.set_zlabel(r"$\log_{10}E$", labelpad=-2)
ax.tick_params(labelsize=7, pad=-1); ax.view_init(elev=38, azim=-125)
ax.set_title("(b) energy of three consecutive zeros", fontsize=10)
fig.subplots_adjust(left=0.07, right=0.97, bottom=0.15, top=0.9, wspace=0.12)
save(fig, "fig_kernel.png")

# Figure 2: the majorant and the eigenvalue bound
d1 = 0.55
fig, axs = plt.subplots(1, 3, figsize=(8, 2.8))
u = np.linspace(-0.9, 0.9, 1000)
ax = axs[0]
ax.plot(u, np.where(np.abs(u) <= 0.5, np.cos(s2*u), 0), color=BLUE, lw=1.6, label=r"$\cos(\sqrt{2}u)\,1_{|u|\leq 1/2}$")
ax.plot(u, np.sinc(d1*u)**2, color=RED, lw=1.6, label=r"$\mathrm{sinc}^2(\pi\delta_1u)$")
ax.set_xlabel("$u$"); ax.set_ylim(0, 1.12); ax.legend(frameon=False, fontsize=7, loc="lower center")
ax.set_title("(a) majorant", fontsize=10)
ax = axs[1]
xi = np.linspace(-1, 1, 1000)
ax.fill_between(xi, np.clip(1 - np.abs(xi)/d1, 0, None)/d1, color=RED, alpha=0.25)
ax.plot(xi, np.clip(1 - np.abs(xi)/d1, 0, None)/d1, color=RED, lw=1.6)
for n in (-1, 1): ax.axvline(n*d1, color=GREY, ls="--", lw=0.8)
ax.set_xlabel(r"$\xi$"); ax.set_title("(b) its transform, support $[-\\delta_1,\\delta_1]$", fontsize=10)
ax = axs[2]
dl = np.linspace(0.4, 1.0, 61); lam = []
for d in dl:
    y = d*np.arange(400); lam.append(np.linalg.eigvalsh(k(y[:, None] - y[None, :]))[-1])
ax.plot(dl, lam, "o", ms=2.5, color=BLUE, label="lattice, 400 points")
ax.plot(dl, 1/(dl*K0), color=RED, lw=1.2, label=r"$1/(\delta K_0)$")
ax.axhline(2, color=GREY, lw=0.8); ax.axvline(d1, color=ORANGE, lw=1.0)
ax.set_xlabel(r"spacing $\delta$"); ax.set_ylabel(r"$\lambda_{\max}$"); ax.legend(frameon=False, fontsize=7)
ax.set_title("(c) largest eigenvalue", fontsize=10)
fig.subplots_adjust(left=0.05, right=0.99, bottom=0.18, top=0.88, wspace=0.32)
save(fig, "fig_majorant.png")

# Figure 3: a block of zeros, its Gram matrix, and the pinching
y = np.cumsum([0, 1.04, 1.98, 1.03, 0.31, 1.07, 1.96, 0.42, 0.29, 1.05, 2.01, 1.04])
n = len(y); G = k(y[:, None] - y[None, :])
blocks = [(0, 3, "run", BLUE), (3, 5, "pair", RED), (5, 6, "run", BLUE), (6, 8, "pair", RED), (8, 9, "single", GREY), (9, 12, "run", BLUE)]
def blk(i):
    for b, (a, c, _, _) in enumerate(blocks):
        if a <= i < c: return b
fig = plt.figure(figsize=(8, 3.5))
ax = fig.add_subplot(1, 2, 1, projection="3d")
order = sorted(((i, j) for i in range(n) for j in range(n)), key=lambda t: t[0] - t[1])
for i, j in order:
    col = blocks[blk(i)][3] if blk(i) == blk(j) else "#cfcfcf"
    ax.bar3d(i, j, 0, 0.75, 0.75, abs(G[i, j]), color=col, shade=True, linewidth=0)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zlim(0, 1); ax.set_zlabel(r"$|k(x_i-x_j)|$", labelpad=-3)
ax.tick_params(labelsize=7, pad=-1); ax.view_init(elev=30, azim=-45)
ax.set_title("(a) Gram matrix; the kept blocks in colour", fontsize=10)
ax = fig.add_subplot(1, 2, 2)
for b, (a, c, lab, col) in enumerate(blocks):
    if c - a > 1: ax.plot([y[a], y[c - 1]], [0, 0], color=col, lw=4, alpha=0.35)
    ax.plot(y[a:c], np.zeros(c - a), "o", ms=7, color=col)
    ax.annotate(lab, ((y[a] + y[c - 1])/2, 0.14 if b % 2 == 0 else -0.2), ha="center", fontsize=9, color=col)
for i in range(n - 1):
    if y[i + 1] - y[i] < 0.55:
        ax.annotate(f"{y[i+1]-y[i]:.2f}", ((y[i] + y[i + 1])/2, 0.06 if i != 6 else -0.09), ha="center", fontsize=7, color=ORANGE)
ax.set_ylim(-0.35, 0.3); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("normalised ordinates of twelve consecutive simple zeros")
ax.set_title("(b) runs, pairs and single zeros", fontsize=10)
fig.subplots_adjust(left=0.02, right=0.98, bottom=0.17, top=0.9, wspace=0.1)
save(fig, "fig_blocks.png")

# Figure 4: what the method gives and where it stops
mu = np.array([0.0015, 0.00175, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.00325, 0.0035, 0.00375, 0.004, 0.00425, 0.0045, 0.00475, 0.005, 0.00525, 0.0055, 0.00575, 0.006, 0.00625, 0.0065, 0.00675, 0.007, 0.00725, 0.0075, 0.00775, 0.008, 0.00825, 0.0085, 0.00875, 0.009])
e = np.array([0.0033696, 0.0038696, 0.0043102, 0.0047403, 0.0051570, 0.0055737, 0.0059725, 0.0063475, 0.0067226, 0.0070976, 0.0074727, 0.0078478, 0.0082228, 0.0085979, 0.0089729, 0.0093480, 0.0097230, 0.0100980, 0.0104731, 0.0108271, 0.0111605, 0.0114939, 0.0118273, 0.0121607, 0.0124941, 0.0128275, 0.0131609, 0.0134943, 0.0138252, 0.0141378, 0.0144503])   # all pairs, best periodic pattern
fig, axs = plt.subplots(1, 2, figsize=(8, 3.1), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axs[0]
names = ["Levinson\n1974", "Conrey\n1989", "Bui, Conrey,\nYoung 2011", "Pratt et al.\n2020", "Alpöge,\nFurman 2026", "this\npaper"]
vals = [1/3, 2/5, 0.41, 5/12, H, 0.6731134]
shown = ["1/3", "2/5", "41%", "5/12", "67.25%", "67.31%"]
cols = [GREY]*4 + [BLUE, RED]
ax.bar(range(len(vals)), vals, color=cols, width=0.65)
for i, v in enumerate(vals): ax.text(i, v + 0.01, shown[i], ha="center", fontsize=7.5)
ax.set_xticks(range(len(vals))); ax.set_xticklabels(names, fontsize=6.5); ax.set_ylim(0, 0.8)
ax.set_ylabel("proportion of zeros"); ax.set_title("(a) zeros on the line: all (grey), simple (colour)", fontsize=10)
ax = axs[1]
ax.plot(mu*1e3, (H - mu)/(1 - e), "s-", color=GREY, ms=2.2, lw=1, label="chain limit of the method")
ax.plot([7/2000*1e3], [0.6730728], "o", color=BLUE, ms=5, label=r"$R=8$, $\mu_0=1/2000$")
ax.plot([7/3200*1e3], [0.6731134], "o", color=RED, ms=5, label=r"$R=8$, $\mu_0=1/3200$")
ax.axhline(0.6730085, color=GREEN, lw=0.9, ls="--"); ax.text(4.3, 0.67303, "seven-point bound", fontsize=7, color=GREEN)
ax.axhline(H, color=BLUE, lw=0.9, ls=":"); ax.text(2.0, 0.67253, r"$H_{\mathrm{MT}}$", fontsize=8, color=BLUE)
ax.set_xlabel(r"$\mu=(R-1)\mu_0$ ($\times10^{-3}$)"); ax.set_ylim(0.6724, 0.6737)
ax.legend(frameon=False, fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.3)); ax.set_title("(b) the windows and the chain", fontsize=10)
ax.ticklabel_format(axis="y", useOffset=False)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.2, top=0.9, wspace=0.32)
save(fig, "fig_limits.png")
