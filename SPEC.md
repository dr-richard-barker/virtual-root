# Virtual Root — model specification (v1)

Extracted from the primary literature to rebuild an open "SimuPlant"/Virtual Root auxin
transport simulator. This is the implementation contract: build to *this*, not to a
guessed geometry.

## Sources
- **Grieneisen et al. 2007**, *Auxin transport is sufficient to generate a maximum and gradient guiding root growth*, [Nature 449:1008](https://www.nature.com/articles/nature06215) — the reflux model SimuPlant is based on.
- **Band et al. 2014**, *Systems Analysis of Auxin Transport in the Arabidopsis Root Apex*, [Plant Cell 26:862](https://academic.oup.com/plcell/article/26/3/862/6099879) — full cell+apoplast model, PIN/AUX1 maps.
- **Review with explicit equations**: *Mathematical Modelling of Auxin Transport in Plant Tissues: Flux Meets Signalling and Growth*, [PMC6976557](https://pmc.ncbi.nlm.nih.gov/articles/PMC6976557/) (Eq. 6).
- **SimuPlant protocol**: *The Virtual Root … using SimuPlant* (2021), [PubMed 34822153](https://pubmed.ncbi.nlm.nih.gov/34822153/).

---

## 1. Compartment model (the key correction)

Auxin lives in **two compartment types**, not one:
1. **Cell cytoplasm** `a_i` — one value per cell.
2. **Apoplast / cell-wall segments** `A_ij` — one value per membrane interface between
   adjacent cells i and j.

Transport is **cell → its wall → neighbour cell**, plus **wall ↔ wall diffusion**. The
2D root tip is a lattice of cells (each cell a box; walls are the interfaces). Per
Grieneisen, diffusion and membrane permeability are **separate** parameters — this is
exactly what my v0 MVP lacked (it had cell-to-cell flux only).

---

## 2. Governing equations

### Rigorous form (Band 2014 / review Eq. 6)

Cytoplasm of cell *i*:
```
da_i/dt = α_a − μ_a·a_i − (1/V_i) Σ_{j~i} S^m_ij · J^a_ij        (+ optional TIR1 terms γ_a·c_i − β_a·a_i·s_i)
```
Apoplast segment between *i* and *j*:
```
dA_ij/dt = (1/V_ij)( S^m_ij·J^a_ij − S^w_ij·J^A_ij − Σ_k S_ijk·J_ijk ) − μ_a·A_ij
```
Trans-membrane flux (cytoplasm i → wall ij), passive + PIN (efflux) + AUX1 (influx),
PIN/AUX1 terms saturating (Michaelis–Menten):
```
J^a_ij = φ_a·κ^ef_a·(a_i − κ^in_a·A_ij)
       + φ_p·P_ij·[ κ^ef_p·a_i/(θ_ap + a_i) − κ^in_p·A_ij/(θ_ap + A_ij) ]
       + φ_u·U_ij·[ κ^ef_u·a_i/(θ_au + a_i) − κ^in_u·A_ij/(θ_au + A_ij) ]
```
- `S^m_ij` membrane interface area; `S^w_ij` wall–wall interface area; `V` volumes.
- `P_ij`, `U_ij` = membrane-bound PIN / AUX1 amount on that side (the localisation maps, §4).
- `J_ijk` = apoplast→apoplast diffusion to neighbouring wall segment k.

### Simpler form (Grieneisen 2007 — recommended for the rebuild)

Per-membrane-**side** permeabilities; no saturation. For each cell side facing wall `w`:
```
efflux  (cell→wall):  p_ef(side) · a_i
influx  (wall→cell):  p_in(side) · A_w
```
```
da_i/dt = prod_i − k·a_i − (1/V_cell) Σ_sides Sm·( p_ef(side)·a_i − p_in(side)·A_w )
dA_w/dt = (1/V_wall)[ Σ_{cells c on w} Sm·( p_ef·a_c − p_in·A_w )
                       + D · Σ_{adjacent walls w'} (A_w' − A_w) ] − k·A_w
```
`p_ef`/`p_in` take the carrier-dependent values in §3; `D` is apoplastic diffusion.

---

## 3. Parameters

### Grieneisen 2007 permeabilities (use these for the rebuild — units µm·s⁻¹)
| Quantity | Value | Notes |
|---|---|---|
| Influx permeability `p_in` | **20** | baseline; **×4 (≈80)** where AUX1 over-expressed/present |
| Efflux permeability `p_ef` | **1 / 5 / 20** | no PIN / weak PIN / strong PIN on that membrane side |
| Apoplastic diffusion `D` | separate term | per Kramer et al. 2007 measurements |
| Production `prod` | baseline + **elevated in QC & columella initials** | the local source that helps seat the maximum |
| Decay `k` | first-order | uniform |

### Review Eq.6 parameter set (units µm·min⁻¹ / µM — alternative, saturating model)
| Symbol | Value | Meaning |
|---|---|---|
| φ_a | 0.55 µm·min⁻¹ | passive membrane permeability |
| φ_p | 0.27 µm·min⁻¹ | PIN saturating flux rate |
| φ_u | 0.55 µm·min⁻¹ | AUX1 saturating flux rate |
| θ_ap = θ_au | 1 µM | half-saturation constants |
| α_a | 0.5 µM·min⁻¹ | auxin biosynthesis |
| μ_a | 0.5 min⁻¹ | auxin degradation |
| κ^ef, κ^in | dimensionless | efflux/influx coefficients |

> ⚠️ **Unit consistency:** Grieneisen uses µm·s⁻¹; the review uses µm·min⁻¹. Pick one
> system and convert everything (×60) before mixing.

---

## 4. Carrier localisation maps (Band 2014)

Guiding principle: **AUX1/LAX (non-polar influx) sets *which* tissues hold auxin; PIN
(polar efflux) sets the *direction* of flow.**

### PIN polarity by tissue
| Membrane orientation | Tissues |
|---|---|
| **Rootward** (basal, ↓) | stele/vasculature, endodermis, meristem cortex, proximal-meristem epidermis |
| **Shootward** (apical, ↑) | lateral root cap, elongation-zone cortex, elongation-zone & distal-meristem epidermis |
| **Omnidirectional** (all faces) | **QC, columella initials, columella tiers S1 & S2** |
| **Lateral, inward** (outer→stele) | outer layers near the tip — closes the reflux loop |

### AUX1 / LAX influx by tissue
| Carrier | Tissues |
|---|---|
| AUX1 | lateral root cap, elongation-zone epidermis, columella S1–S3, elongating cortex |
| LAX2 | QC, columella initials, rootward half of meristematic stele |
| LAX3 | columella tier S2 only |

### Production
Baseline in all cells; **higher in the QC and columella initials**.

---

## 5. Why my v0 MVP failed → the fixes
| v0 problem | Fix from this spec |
|---|---|
| Single compartment (cell-to-cell flux) | Add apoplast wall compartments; flux = cell↔wall + wall↔wall diffusion (§1) |
| One combined transport coefficient | Separate `p_in`, `p_ef` per side + diffusion `D` (§2–3) |
| Guessed rectangular PIN map | Use the Band 2014 PIN map incl. **omnidirectional QC/columella** + **lateral-inward** outer layers (§4) |
| No AUX1 trap that actually holds auxin | AUX1/LAX influx in LRC + columella (high `p_in`) traps auxin at the tip (§4) |
| No local source at the tip | **Elevated production in QC + columella** (§3–4) |
| Auxin piled at boundaries | Closed reflux loop (lateral-inward return) + the AUX1 trap → maximum seats at QC |

## 6. Validation target
Steady state must show the **auxin maximum at the QC / columella initials**, a gradient
decreasing shootward, and visible reflux (down the stele, around the tip, up the outer
layers). Reproduce Grieneisen 2007 Fig. 1–2 qualitatively before adding the UI.

## 7. Build order
1. Implement the §2 *Grieneisen* form on a small **digitised root cell template** (not a
   plain rectangle) with the §4 maps and §3 parameters.
2. Confirm the §6 validation target.
3. Add interactive controls (PIN/AUX1 strength, `D`, decay, production) + PNG/data export.
4. Ship as pure-JS on GitHub Pages (durable) and embed in AIRI Stage VII.

## 8. Validating the model against DII-VENUS (for researchers)

The model predicts a per-cell auxin distribution; the **DII-VENUS** biosensor lets you
test that prediction against a living root.

**What it is.** DII-VENUS ([Brunoud et al. 2012, *Nature*](https://doi.org/10.1038/nature10791))
fuses the auxin-binding domain II (DII) of an Aux/IAA to fast-maturing VENUS. Auxin drives
its TIR1/AFB-dependent degradation, so **fluorescence is *inversely* related to auxin** —
bright cells = low auxin, dark cells = high auxin. It reports the auxin *input* directly and
fast (unlike the slower transcriptional DR5 *output* reporter).

**Validation workflow**
1. **Image** the root tip by confocal microscopy. Prefer the ratiometric **R2D2**
   ([Liao et al. 2015, *Nat. Methods*](https://doi.org/10.1038/nmeth.3279)), which adds a
   non-degradable *mDII-tdTomato* control so an mDII/DII ratio corrects for expression level
   and imaging depth.
2. **Segment** cells and measure per-cell mean fluorescence.
3. **Convert fluorescence → relative auxin.** Because the signal is inverse and non-linear,
   invert it using the DII degradation model (auxin-dependent degradation via TIR1 binding),
   as in [Band et al. 2012](https://doi.org/10.1073/pnas.1201498109) /
   [2014](https://doi.org/10.1105/tpc.113.119495). R2D2 ratios give a more direct readout.
4. **Compare cell-by-cell** with the model's steady state:
   - Does the **maximum** sit at the QC / columella initials? *(the tool reproduces this)*
   - Are the **shootward gradient** and the **lateral-root-cap** signal reproduced?
   - Fit model parameters (permeabilities, production) to the measured pattern.
5. **Dynamic test — gravitropism.** Gravistimulate and time-lapse DII-VENUS/R2D2: a lateral
   auxin asymmetry develops on the lower side within minutes
   ([Band et al. 2012](https://doi.org/10.1073/pnas.1201498109)). The tool's
   **Gravistimulate** preset predicts exactly this — compare the predicted vs measured
   lower/upper ratio *and its time course*.

**Caveats.** DII-VENUS reports auxin *perception*, not transcriptional output; the signal is
inverse and saturates; VENUS maturation time and photobleaching bias fast dynamics. Use R2D2
for quantitative work and calibrate the fluorescence→auxin conversion for each imaging setup.

**References:** Brunoud *et al.* (2012) *Nature* 482:103; Band *et al.* (2012) *PNAS* 109:4668;
Band *et al.* (2014) *Plant Cell* 26:862; Liao *et al.* (2015) *Nature Methods* 12:207.
