"""
Virtual Root v2 — model-fidelity polish over v1:
  (1) rounded root-tip template (domed cap, not a rectangle),
  (2) smoother shootward gradient (a shoot-derived auxin supply fills the root, so auxin
      is present throughout with the QC still the maximum), and
  (3) physical units: Grieneisen permeabilities in um/s, cell size in um, time in seconds,
      concentrations in uM.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------- geometry (rounded tip) ----------
R, C = 52, 17
c0 = C // 2
r_qc, r_mer, CAP = 46, 30, 4.7
rows = np.arange(R)[:, None].repeat(C, 1)
rad = np.abs(np.arange(C) - c0)[None, :].repeat(R, 0)
active = np.zeros((R, C), bool)
for r in range(R):
    for c in range(C):
        rd = abs(c - c0)
        if r < r_qc:
            active[r, c] = rd <= 4
        else:                                    # domed cap: quarter-circle taper
            active[r, c] = rd <= np.sqrt(max(0.0, CAP**2 - (r - r_qc)**2))

stele  = rad <= 1; endo = rad == 2; cortex = rad == 3; epi = rad == 4
tip    = rows >= r_qc
qc        = active & (rows == r_qc) & (rad <= 1)
columella = active & (rows > r_qc) & (rad <= 2)
lrc       = active & (rows >= r_qc - 2) & (rad >= 3)

# ---------- physical parameters ----------
L = 10.0                      # cell size (um)
g = 1.0 / L                   # surface/volume factor (1/um) in 2D
EF_BASE, EF_PIN = 1.0, 16.0   # efflux permeability (um/s): none / strong PIN
IN_BASE, IN_AUX1 = 20.0, 42.0 # influx permeability (um/s): baseline / AUX1 (moderate trap)
dt = 0.04                     # s
k  = 0.0012                   # auxin decay (1/s)  -> lets a gradient form
# Production (uM/s): a strong shoot supply fills the root; the tip is the maximum because
# AUX1 traps auxin there, with a small local QC/columella source. A modest baseline keeps
# auxin present throughout (a visible shootward gradient).
S_shoot, P_qc, P_col, P_base = 0.30, 0.025, 0.012, 0.005

# ---------- carrier maps ----------
ef_up = np.full((R, C), EF_BASE); ef_dn = np.full((R, C), EF_BASE)
ef_lf = np.full((R, C), EF_BASE); ef_rt = np.full((R, C), EF_BASE)
rootward  = ((stele | endo) & ~tip) | (cortex & (rows >= r_mer) & ~tip)
shootward = (cortex & (rows < r_mer)) | (epi & (rows < r_mer)) | (epi & tip) | lrc
ef_dn[rootward] = EF_PIN
ef_up[shootward] = EF_PIN
inward = epi & (rows >= r_mer)
ef_rt[inward & (np.arange(C)[None, :] < c0)] = EF_PIN
ef_lf[inward & (np.arange(C)[None, :] > c0)] = EF_PIN
for arr in (ef_up, ef_dn, ef_lf, ef_rt): arr[columella] = EF_PIN   # columella redistributes

pin_in = np.full((R, C), IN_BASE)
pin_in[epi | lrc | columella | qc | (cortex & tip)] = IN_AUX1

prod = np.where(active, P_base, 0.0)
prod[columella] = P_col
prod[qc] = P_qc
prod[0, c0-1:c0+2] = S_shoot              # shoot-derived auxin enters the top of the stele

for arr in (ef_up, ef_dn, ef_lf, ef_rt, pin_in, prod): arr[~active] = 0.0

# ---------- integrate ----------
a = np.where(active, 0.1, 0.0)
for step in range(200000):
    bV = active[:, :-1] & active[:, 1:]; dV = pin_in[:, :-1] + pin_in[:, 1:]
    Wv = np.where(bV, (ef_rt[:, :-1]*a[:, :-1] + ef_lf[:, 1:]*a[:, 1:]) / np.where(dV>0, dV, 1), 0)
    fL = (pin_in[:, :-1]*Wv - ef_rt[:, :-1]*a[:, :-1]) * bV
    fR = (pin_in[:, 1:]*Wv - ef_lf[:, 1:]*a[:, 1:]) * bV
    bH = active[:-1, :] & active[1:, :]; dH = pin_in[:-1, :] + pin_in[1:, :]
    Wh = np.where(bH, (ef_dn[:-1, :]*a[:-1, :] + ef_up[1:, :]*a[1:, :]) / np.where(dH>0, dH, 1), 0)
    fU = (pin_in[:-1, :]*Wh - ef_dn[:-1, :]*a[:-1, :]) * bH
    fD = (pin_in[1:, :]*Wh - ef_up[1:, :]*a[1:, :]) * bH
    da = prod - k*a
    da[:, :-1] += g*fL; da[:, 1:] += g*fR; da[:-1, :] += g*fU; da[1:, :] += g*fD
    # open-boundary efflux (auxin pumped toward soil/shoot leaves the domain)
    aL = np.zeros_like(active); aL[:,1:]=active[:,:-1]; aR=np.zeros_like(active); aR[:,:-1]=active[:,1:]
    aU = np.zeros_like(active); aU[1:,:]=active[:-1,:]; aD=np.zeros_like(active); aD[:-1,:]=active[1:,:]
    da[active & ~aL] -= g*ef_lf[active & ~aL]*a[active & ~aL]
    da[active & ~aR] -= g*ef_rt[active & ~aR]*a[active & ~aR]
    da[active & ~aU] -= g*ef_up[active & ~aU]*a[active & ~aU]
    da[active & ~aD] -= g*ef_dn[active & ~aD]*a[active & ~aD]
    a += dt*da; a[~active] = 0.0
    if step % 1000 == 0 and np.max(np.abs(da[active])) < 1e-8: break

am = np.where(active, a, -1); mr, mc = np.unravel_index(np.argmax(am), am.shape)
print(f"step {step}, sim-time {step*dt:.0f}s, max auxin {a[mr,mc]:.3f} uM at (row={mr}, col={mc})")
print(f"QC row {r_qc}, axis col {c0}; max at QC: {mr>=r_qc-1 and abs(mc-c0)<=2}")
# how far up does the gradient reach (fraction of stele cells above the meristem with >5% of max)
upper = a[:r_mer, c0] / a[mr,mc]
print(f"auxin in upper stele (frac of max): top={upper[0]:.2f}, mid={upper[r_mer//2]:.2f}")

from matplotlib.colors import PowerNorm
disp = np.where(active, a, np.nan)
fig, ax = plt.subplots(figsize=(5.6, 8.2))
im = ax.imshow(disp, cmap="magma", origin="upper", aspect="equal",
               norm=PowerNorm(gamma=0.45, vmin=0, vmax=np.nanmax(disp)))
ax.plot(mc, mr, "o", mfc="none", mec="lime", ms=13, mew=2)

# --- flux arrows (efflux polarity x concentration), subsampled ---
vx = (ef_rt - ef_lf) * a; vy = (ef_dn - ef_up) * a
xs, ys, us, vs = [], [], [], []
for r in range(2, R, 3):
    for c in range(1, C, 2):
        if active[r, c] and (abs(vx[r, c]) + abs(vy[r, c])) > 0.5:
            mg = np.hypot(vx[r, c], vy[r, c]); xs.append(c); ys.append(r)
            us.append(vx[r, c]/mg); vs.append(vy[r, c]/mg)
ax.quiver(xs, ys, us, vs, color="white", alpha=0.55, scale=26, width=0.005, headwidth=4, pivot="mid")

# --- anatomy labels ---
L_ = dict(fontsize=8, color="black", arrowprops=dict(arrowstyle="-", color="0.45", lw=0.6))
ax.annotate("Elongation\nzone", xy=(c0-4, 14), xytext=(-10, 14), va="center", ha="right", **L_)
ax.annotate("Meristem", xy=(c0-4, 32), xytext=(-10, 32), va="center", ha="right", **L_)
ax.annotate("QC", xy=(mc, mr), xytext=(-10, mr), va="center", ha="right", fontsize=8, color="lime",
            arrowprops=dict(arrowstyle="-", color="lime", lw=0.6))
ax.annotate("Columella\n/ root cap", xy=(c0, r_qc+2), xytext=(-10, r_qc+3), va="center", ha="right", **L_)
ax.annotate("Stele", xy=(c0, 20), xytext=(C+5, 20), va="center", **L_)
ax.annotate("Cortex", xy=(c0+3, 24), xytext=(C+5, 24), va="center", **L_)
ax.annotate("Epidermis", xy=(c0+4, 28), xytext=(C+5, 28), va="center", **L_)

# --- demo bend arrow at the tip (placement check for the JS version) ---
ax.annotate("", xy=(c0-3.5, R+1.5), xytext=(c0, R+1.5), arrowprops=dict(arrowstyle="->", color="#0088bb", lw=2.2))
ax.text(c0-1.7, R+3, "bends toward\nlower side", color="#0088bb", fontsize=7, ha="center", va="top")

ax.set_xlim(-12, C+9); ax.set_ylim(R+5, -2)
ax.set_title("Virtual Root — anatomy, auxin & reflux", fontsize=10)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("← radial →"); ax.set_ylabel("shoot → tip")
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="auxin (uM)")
fig.tight_layout(); fig.savefig(__file__.replace("auxin_v2.py","auxin_v2_annotated.png"), dpi=130)
fig.savefig(__file__.replace("auxin_v2.py","auxin_v2.png"), dpi=130)
print("saved auxin_v2_annotated.png")
