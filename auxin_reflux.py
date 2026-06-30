"""
Virtual Root (open SimuPlant) — MVP kernel.

A cell-grid model of polar auxin transport in the Arabidopsis root tip, after the
Grieneisen et al. (2007) reflux / reverse-fountain model. Each cell exchanges auxin
with its 4 neighbours by:
  - passive diffusion across the membrane (permeability D), and
  - active PIN-mediated efflux (rate p) in a tissue-specific polar direction.
Auxin enters from the shoot end of the stele (source S) and is removed by decay (k).

Goal of this MVP: show that prescribing PIN polarity per tissue layer is sufficient to
self-organise an auxin maximum at the quiescent centre (QC) — the hallmark result.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---- geometry ----
R, C = 44, 21            # rows (0=shoot end, R-1=tip), cols
c0 = C // 2              # central axis
r_tip = 36               # rows >= r_tip are the tip zone (columella / root cap / QC)
STELE = 1                # stele half-width: |c-c0| <= STELE

# ---- parameters ----
D = 0.20                 # passive permeability (diffusion)
p = 0.70                 # PIN active-transport (efflux) strength
a = 1.20                 # AUX1 active-influx strength (tip trap)
k = 0.02                 # auxin decay rate
S = 0.50                 # auxin influx at the shoot-end source cells
dt = 0.10
STEPS, TOL = 40000, 1e-7

# ---- PIN polarity field: fraction of each cell's PIN on up/down/left/right membrane ----
pin_up = np.zeros((R, C)); pin_dn = np.zeros((R, C))
pin_lf = np.zeros((R, C)); pin_rt = np.zeros((R, C))

for r in range(R):
    for c in range(C):
        in_stele = abs(c - c0) <= STELE
        if r < r_tip:
            if in_stele:
                pin_dn[r, c] = 1.0          # stele pumps auxin rootward (down)
            else:
                pin_up[r, c] = 1.0          # outer layers pump auxin shootward (up & out)
        else:                                # tip zone (columella / root cap)
            if in_stele:                     # columella below QC: reflect auxin UP toward QC,
                pin_up[r, c] = 0.8           # leaking a little laterally to the cap
                pin_lf[r, c] = pin_rt[r, c] = 0.1
            else:                            # lateral root cap: return auxin shootward
                pin_up[r, c] = 1.0

# ---- AUX1 influx field: high in the tip (columella + root cap) to trap auxin there ----
aux1 = np.zeros((R, C))
aux1[r_tip:, :] = 1.0                 # tip zone traps auxin (the QC auxin maximum)
aux1[r_tip - 3:r_tip, :] = 0.5        # graded just above the tip

# ---- source (shoot-derived auxin enters top of stele) ----
prod = np.zeros((R, C))
prod[0, c0 - STELE:c0 + STELE + 1] = S

# ---- integrate to steady state ----
A = np.zeros((R, C))
for step in range(STEPS):
    # vertical net flow from (r,c) -> (r+1,c)
    Tv = (D * (A[:-1, :] - A[1:, :])
          + p * pin_dn[:-1, :] * A[:-1, :]
          - p * pin_up[1:, :] * A[1:, :]
          + a * aux1[1:, :] * A[:-1, :]          # AUX1 in lower cell pulls auxin down
          - a * aux1[:-1, :] * A[1:, :])         # AUX1 in upper cell pulls auxin up
    # horizontal net flow from (r,c) -> (r,c+1)
    Th = (D * (A[:, :-1] - A[:, 1:])
          + p * pin_rt[:, :-1] * A[:, :-1]
          - p * pin_lf[:, 1:] * A[:, 1:]
          + a * aux1[:, 1:] * A[:, :-1]          # AUX1 in right cell pulls auxin right
          - a * aux1[:, :-1] * A[:, 1:])         # AUX1 in left cell pulls auxin left
    dA = prod - k * A
    dA[:-1, :] -= Tv; dA[1:, :] += Tv
    dA[:, :-1] -= Th; dA[:, 1:] += Th
    # open-boundary active efflux: auxin pumped toward a domain edge leaves the system
    dA[0, :]  -= p * pin_up[0, :]  * A[0, :]
    dA[-1, :] -= p * pin_dn[-1, :] * A[-1, :]
    dA[:, 0]  -= p * pin_lf[:, 0]  * A[:, 0]
    dA[:, -1] -= p * pin_rt[:, -1] * A[:, -1]
    A += dt * dA
    if step % 200 == 0 and np.max(np.abs(dA)) < TOL:
        break

# ---- report ----
qc_r, qc_c = np.unravel_index(np.argmax(A), A.shape)
print(f"converged at step {step}")
print(f"auxin max = {A.max():.3f} at (row={qc_r}, col={qc_c}); axis col={c0}, tip starts row {r_tip}")
print(f"max is in the tip zone: {qc_r >= r_tip - 1}; on the central axis: {abs(qc_c - c0) <= STELE}")

# ---- visualise ----
fig, ax = plt.subplots(figsize=(4.2, 7))
im = ax.imshow(A, cmap="magma", origin="upper", aspect="equal")
ax.axvspan(c0 - STELE - 0.5, c0 + STELE + 0.5, color="cyan", alpha=0.08)
ax.axhline(r_tip - 0.5, color="white", lw=0.6, ls="--", alpha=0.5)
ax.plot(qc_c, qc_r, "o", mfc="none", mec="lime", ms=14, mew=2)
ax.annotate("auxin maximum", (qc_c, qc_r), (qc_c + 2, qc_r),
            color="lime", fontsize=8, va="center")
ax.set_title("Virtual Root: auxin reflux model\n(stele↓  outer↑  tip redistributes)", fontsize=9)
ax.set_xlabel("← radial →"); ax.set_ylabel("shoot end  →  root tip")
ax.set_xticks([]); ax.set_yticks([])
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="auxin concentration")
fig.tight_layout()
out = __file__.replace("auxin_reflux.py", "auxin_reflux.png")
fig.savefig(out, dpi=130)
print("saved", out)
