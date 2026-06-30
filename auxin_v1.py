"""
Virtual Root v1 — two-compartment (cell + apoplast wall) auxin transport, built to SPEC.md.

Cells (a) and the wall segments between them (W) are separate compartments. Each cell
side has an EFFLUX permeability (low=1 baseline, high=20 where PIN points that way) and
each cell has an INFLUX permeability (20 baseline, 80 where AUX1 present) — the Grieneisen
(2007) scheme. Walls are treated at quasi-steady state. Tissue-specific PIN polarity +
AUX1 influx + elevated QC/columella production should seat the auxin maximum at the QC.
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------- geometry ----------
R, C = 50, 15
c0 = C // 2
rows = np.arange(R)[:, None].repeat(C, 1)
rad = np.abs(np.arange(C) - c0)[None, :].repeat(R, 0)     # radial index from axis
active = rad <= 4                                          # root = stele..epidermis (width 9)

r_qc, r_mer = 44, 30                                       # QC row; meristem upper bound
stele  = rad <= 1
endo   = rad == 2
cortex = rad == 3
epi    = rad == 4
tip    = rows >= r_qc
qc        = active & (rows == r_qc) & (rad <= 1)
columella = active & (rows > r_qc) & (rad <= 2)
lrc       = active & (rows >= r_qc - 2) & (rad >= 3)        # lateral root cap (outer tip)
omni      = qc | columella

# ---------- carrier maps (SPEC §3-4) ----------
EF_BASE, EF_PIN, IN_BASE, IN_AUX1 = 1.0, 20.0, 20.0, 80.0
ef_up = np.full((R, C), EF_BASE); ef_dn = np.full((R, C), EF_BASE)
ef_lf = np.full((R, C), EF_BASE); ef_rt = np.full((R, C), EF_BASE)

# PIN rootward (down): stele, endodermis, meristem cortex (above the tip)
rootward = ((stele | endo) & ~tip) | (cortex & (rows >= r_mer) & ~tip)
ef_dn[rootward] = EF_PIN
# PIN shootward (up): elongation cortex + epidermis (above meristem), and the cap returns up
shootward = (cortex & (rows < r_mer)) | (epi & (rows < r_mer)) | (epi & tip) | lrc
ef_up[shootward] = EF_PIN
# PIN lateral-inward in the tip epidermis -> closes the reflux loop back to the stele
inward_tip = epi & (rows >= r_mer)
ef_rt[inward_tip & (np.arange(C)[None, :] < c0)] = EF_PIN
ef_lf[inward_tip & (np.arange(C)[None, :] > c0)] = EF_PIN
# PIN omnidirectional in the COLUMELLA only (it redistributes auxin laterally).
# The QC itself stays low-efflux -> a trap where downward stele flux accumulates.
for arr in (ef_up, ef_dn, ef_lf, ef_rt):
    arr[columella] = EF_PIN

# AUX1 influx (non-polar): epidermis, cap, columella, QC, tip cortex
pin_in = np.full((R, C), IN_BASE)
aux1 = epi | lrc | columella | qc | (cortex & tip)
pin_in[aux1] = IN_AUX1

# production: baseline everywhere, elevated in columella, highest at the QC
prod = np.where(active, 0.1, 0.0)
prod[columella] = 1.0
prod[qc] = 2.0

k = 0.05                                                   # decay
for arr in (ef_up, ef_dn, ef_lf, ef_rt, pin_in, prod):
    arr[~active] = 0.0

# neighbour-active masks (for boundary efflux to "shoot"/soil)
actL = np.zeros_like(active); actL[:, 1:] = active[:, :-1]
actR = np.zeros_like(active); actR[:, :-1] = active[:, 1:]
actU = np.zeros_like(active); actU[1:, :] = active[:-1, :]
actD = np.zeros_like(active); actD[:-1, :] = active[1:, :]
bL, bR = active & ~actL, active & ~actR
bU, bD = active & ~actU, active & ~actD

# ---------- integrate (quasi-steady-state walls) ----------
a = np.where(active, 0.2, 0.0)
dt, STEPS, TOL = 0.004, 120000, 1e-7
for step in range(STEPS):
    bothV = active[:, :-1] & active[:, 1:]
    denV = pin_in[:, :-1] + pin_in[:, 1:]
    Wv = np.where(bothV, (ef_rt[:, :-1] * a[:, :-1] + ef_lf[:, 1:] * a[:, 1:]) / np.where(denV > 0, denV, 1), 0)
    fL = (pin_in[:, :-1] * Wv - ef_rt[:, :-1] * a[:, :-1]) * bothV
    fR = (pin_in[:, 1:] * Wv - ef_lf[:, 1:] * a[:, 1:]) * bothV

    bothH = active[:-1, :] & active[1:, :]
    denH = pin_in[:-1, :] + pin_in[1:, :]
    Wh = np.where(bothH, (ef_dn[:-1, :] * a[:-1, :] + ef_up[1:, :] * a[1:, :]) / np.where(denH > 0, denH, 1), 0)
    fU = (pin_in[:-1, :] * Wh - ef_dn[:-1, :] * a[:-1, :]) * bothH
    fD = (pin_in[1:, :] * Wh - ef_up[1:, :] * a[1:, :]) * bothH

    da = prod - k * a
    da[:, :-1] += fL; da[:, 1:] += fR
    da[:-1, :] += fU; da[1:, :] += fD
    # boundary efflux: auxin pumped toward an edge/soil leaves the system
    da[bL] -= ef_lf[bL] * a[bL]; da[bR] -= ef_rt[bR] * a[bR]
    da[bU] -= ef_up[bU] * a[bU]; da[bD] -= ef_dn[bD] * a[bD]

    a += dt * da
    a[~active] = 0.0
    if step % 500 == 0 and np.max(np.abs(da[active])) < TOL:
        break

# ---------- report ----------
am = np.where(active, a, -1)
mr, mc = np.unravel_index(np.argmax(am), am.shape)
print(f"converged at step {step}, max auxin {a[mr,mc]:.3f} at (row={mr}, col={mc})")
print(f"QC at row {r_qc}, axis col {c0}; max in tip zone: {mr >= r_qc-1}; near axis: {abs(mc-c0)<=2}")

# ---------- plot ----------
disp = np.where(active, a, np.nan)
fig, ax = plt.subplots(figsize=(3.6, 7.2))
im = ax.imshow(disp, cmap="magma", origin="upper", aspect="equal")
ax.plot(mc, mr, "o", mfc="none", mec="lime", ms=14, mew=2)
ax.annotate("auxin max", (mc, mr), (mc + 1.5, mr), color="lime", fontsize=8, va="center")
ax.axhline(r_qc - 0.5, color="cyan", lw=0.6, ls="--", alpha=0.6)
ax.set_title("Virtual Root v1 — two-compartment\ncell+wall, PIN/AUX1 maps", fontsize=9)
ax.set_xlabel("← radial →"); ax.set_ylabel("shoot end  →  root tip")
ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04, label="auxin (cytoplasm)")
fig.tight_layout()
out = __file__.replace("auxin_v1.py", "auxin_v1.png")
fig.savefig(out, dpi=130)
print("saved", out)
