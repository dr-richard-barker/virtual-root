"""
Gravity series for the Virtual Root auxin model — µG, Moon, Mars, Earth, 2 g.

The reverse-fountain redistribution needs a *directional* gravity cue: statoliths
sediment, PIN3 (columella) re-localises, and auxin is biased to the lower side. Here the
columella lateral-PIN bias scales with the gravity level g (saturating), so the model
predicts how the lateral auxin asymmetry — the first step of bending — changes from
Earth to the Moon, Mars, microgravity and hypergravity.
Outputs: gravity_series.png (auxin maps) + asymmetry_curve.png (asymmetry vs g).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

# ---- geometry (v2 rounded tip) ----
R, C, c0 = 52, 17, 8
r_qc, r_mer, CAP = 46, 30, 4.7
rows = np.arange(R)[:, None].repeat(C, 1)
rad = np.abs(np.arange(C) - c0)[None, :].repeat(R, 0)
active = np.zeros((R, C), bool)
for r in range(R):
    for c in range(C):
        rd = abs(c - c0)
        active[r, c] = rd <= 4 if r < r_qc else rd <= np.sqrt(max(0.0, CAP**2 - (r - r_qc)**2))
stele = rad <= 1; endo = rad == 2; cortex = rad == 3; epi = rad == 4; tip = rows >= r_qc
qc = active & (rows == r_qc) & (rad <= 1)
col = active & (rows > r_qc) & (rad <= 2)
lrc = active & (rows >= r_qc - 2) & (rad >= 3)

# ---- parameters ----
L = 10.0; g = 1.0 / L; dt = 0.04; k = 0.0012
EF_BASE, EF_PIN, IN_BASE, IN_AUX1 = 1.0, 16.0, 20.0, 42.0
S_shoot, P_qc, P_col, P_base = 0.30, 0.025, 0.012, 0.005
colcols = np.arange(C)[None, :]

def bias(gl):                       # statolith-sedimentation-like saturating response
    return min(0.9, 0.6 * np.tanh(gl) / np.tanh(1.0))

def build(gl, pin37=1.0):
    ef_up = np.full((R, C), EF_BASE); ef_dn = ef_up.copy(); ef_lf = ef_up.copy(); ef_rt = ef_up.copy()
    rootward = ((stele | endo) & ~tip) | (cortex & (rows >= r_mer) & ~tip)
    shootward = (cortex & (rows < r_mer)) | (epi & (rows < r_mer)) | (epi & tip) | lrc
    ef_dn[rootward] = EF_PIN; ef_up[shootward] = EF_PIN
    inw = epi & (rows >= r_mer)
    ef_rt[inw & (colcols < c0)] = EF_PIN; ef_lf[inw & (colcols > c0)] = EF_PIN
    cp = EF_PIN * pin37; b = bias(gl)
    for arr in (ef_up, ef_dn, ef_lf, ef_rt): arr[col] = cp
    ef_lf[col] = cp * (1 + b)        # bias auxin to the LOWER (left) side
    ef_rt[col] = cp * (1 - b)
    pin_in = np.full((R, C), IN_BASE)
    pin_in[epi | lrc | col | qc | (cortex & tip)] = IN_AUX1
    prod = np.where(active, P_base, 0.0); prod[col] = P_col; prod[qc] = P_qc
    prod[0, c0 - 1:c0 + 2] += S_shoot
    for arr in (ef_up, ef_dn, ef_lf, ef_rt, pin_in, prod): arr[~active] = 0.0
    return ef_up, ef_dn, ef_lf, ef_rt, pin_in, prod

def run(gl, pin37=1.0, steps=45000):
    ef_up, ef_dn, ef_lf, ef_rt, pin_in, prod = build(gl, pin37)
    a = np.where(active, 0.1, 0.0)
    aL = np.zeros_like(active); aR = aL.copy(); aU = aL.copy(); aD = aL.copy()
    aL[:, 1:] = active[:, :-1]; aR[:, :-1] = active[:, 1:]; aU[1:, :] = active[:-1, :]; aD[:-1, :] = active[1:, :]
    for _ in range(steps):
        bV = active[:, :-1] & active[:, 1:]; dV = pin_in[:, :-1] + pin_in[:, 1:]
        Wv = np.where(bV, (ef_rt[:, :-1]*a[:, :-1] + ef_lf[:, 1:]*a[:, 1:]) / np.where(dV > 0, dV, 1), 0)
        fL = (pin_in[:, :-1]*Wv - ef_rt[:, :-1]*a[:, :-1]) * bV
        fR = (pin_in[:, 1:]*Wv - ef_lf[:, 1:]*a[:, 1:]) * bV
        bH = active[:-1, :] & active[1:, :]; dH = pin_in[:-1, :] + pin_in[1:, :]
        Wh = np.where(bH, (ef_dn[:-1, :]*a[:-1, :] + ef_up[1:, :]*a[1:, :]) / np.where(dH > 0, dH, 1), 0)
        fU = (pin_in[:-1, :]*Wh - ef_dn[:-1, :]*a[:-1, :]) * bH
        fD = (pin_in[1:, :]*Wh - ef_up[1:, :]*a[1:, :]) * bH
        da = prod - k*a
        da[:, :-1] += g*fL; da[:, 1:] += g*fR; da[:-1, :] += g*fU; da[1:, :] += g*fD
        da[active & ~aL] -= g*ef_lf[active & ~aL]*a[active & ~aL]
        da[active & ~aR] -= g*ef_rt[active & ~aR]*a[active & ~aR]
        da[active & ~aU] -= g*ef_up[active & ~aU]*a[active & ~aU]
        da[active & ~aD] -= g*ef_dn[active & ~aD]*a[active & ~aD]
        a += dt*da; a[~active] = 0.0
    return a

def asym(a):                        # lower(left) - upper(right) auxin in the elongation zone, %
    lo = a[r_mer:r_qc-2, :c0][active[r_mer:r_qc-2, :c0]].sum()
    hi = a[r_mer:r_qc-2, c0+1:][active[r_mer:r_qc-2, c0+1:]].sum()
    return 100 * (lo - hi) / (lo + hi + 1e-9)

LEVELS = [("Microgravity\n(µG, 0 g)", 0.0), ("Moon\n(0.16 g)", 0.16),
          ("Mars\n(0.38 g)", 0.38), ("Earth\n(1 g)", 1.0), ("Hyper-g\n(2 g)", 2.0)]

# ---- Figure 1: auxin maps across gravity ----
fig, axes = plt.subplots(1, len(LEVELS), figsize=(11, 5.4))
results = []
for ax, (name, gl) in zip(axes, LEVELS):
    a = run(gl); results.append((name, gl, a, asym(a)))
    disp = np.where(active, a, np.nan)
    ax.imshow(disp, cmap="magma", origin="upper", aspect="equal",
              norm=PowerNorm(0.45, vmin=0, vmax=np.nanmax(disp)))
    ax.set_title(f"{name}\nasymmetry {asym(a):+.0f}%", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
axes[0].set_ylabel("shoot → tip", fontsize=9)
fig.suptitle("Virtual Root: gravity sets the reverse-fountain auxin asymmetry", fontsize=12)
fig.tight_layout(); fig.savefig("gravity_series.png", dpi=140); print("saved gravity_series.png")

# ---- Figure 2: asymmetry vs gravity ----
gg = np.linspace(0, 2, 13); yy = [asym(run(x, steps=25000)) for x in gg]
fig2, ax2 = plt.subplots(figsize=(6.4, 4.2))
ax2.plot(gg, yy, color="#2e7d32", lw=2)
for name, gl, a, s in results:
    ax2.plot(gl, s, "o", ms=8, color="#d1495b")
    ax2.annotate(name.replace("\n", " "), (gl, s), (gl, s + 3), fontsize=7, ha="center")
ax2.set_xlabel("gravity level (g)"); ax2.set_ylabel("lower-side auxin asymmetry (%)")
ax2.set_title("Predicted gravitropic auxin asymmetry vs gravity")
ax2.axhline(0, color="0.7", lw=0.8); ax2.grid(alpha=0.25)
fig2.tight_layout(); fig2.savefig("asymmetry_curve.png", dpi=140); print("saved asymmetry_curve.png")
for name, gl, a, s in results: print(f"{name.splitlines()[0]:14} g={gl:<4} asymmetry={s:+.1f}%")
