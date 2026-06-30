// Virtual Root — in-browser auxin transport model (port of auxin_v2.py).
// Two compartments (cell + apoplast wall, walls at quasi-steady state); per-side,
// carrier-modulated permeabilities (Grieneisen um/s); Band PIN/AUX1 maps; rounded root
// tip; shoot-derived supply + small QC source -> auxin maximum at the QC atop a gradient.

const R = 52, C = 17, c0 = 8, N = R * C;
const r_qc = 46, r_mer = 30, CAP = 4.7;
const L = 10.0, g = 1.0 / L, dt = 0.04;        // cell size (um), S/V factor, timestep (s)
const EF_BASE = 1, IN_BASE = 20;
const S_shoot = 0.30, P_col = 0.012, P_base = 0.005;   // production (uM/s)

// ---- static tissue masks ----
const active = new Uint8Array(N);
const m = { stele:U(), endo:U(), cortex:U(), epi:U(), qc:U(), col:U(), lrc:U(), tip:U() };
function U(){ return new Uint8Array(N); }
for (let r = 0; r < R; r++) for (let c = 0; c < C; c++) {
  const i = r*C+c, rad = Math.abs(c - c0);
  active[i] = (r < r_qc) ? (rad <= 4 ? 1 : 0)
                         : (rad <= Math.sqrt(Math.max(0, CAP*CAP - (r-r_qc)*(r-r_qc))) ? 1 : 0);
  if (!active[i]) continue;
  m.stele[i]  = rad <= 1 ? 1 : 0;
  m.endo[i]   = rad === 2 ? 1 : 0;
  m.cortex[i] = rad === 3 ? 1 : 0;
  m.epi[i]    = rad === 4 ? 1 : 0;
  m.tip[i]    = r >= r_qc ? 1 : 0;
  m.qc[i]     = (r === r_qc && rad <= 1) ? 1 : 0;
  m.col[i]    = (r > r_qc && rad <= 2) ? 1 : 0;
  m.lrc[i]    = (r >= r_qc - 2 && rad >= 3) ? 1 : 0;
}

// ---- numeric parameter fields ----
const ef_up = new Float64Array(N), ef_dn = new Float64Array(N);
const ef_lf = new Float64Array(N), ef_rt = new Float64Array(N);
const pin_in = new Float64Array(N), prod = new Float64Array(N);
let kDecay = 0.0012;
const params = { EF_PIN: 16, IN_AUX1: 42, P_qc: 0.025, k: 0.0012, gravity: false };

function rebuild() {
  kDecay = params.k;
  for (let r = 0; r < R; r++) for (let c = 0; c < C; c++) {
    const i = r*C+c;
    if (!active[i]) { ef_up[i]=ef_dn[i]=ef_lf[i]=ef_rt[i]=pin_in[i]=prod[i]=0; continue; }
    const rad = Math.abs(c - c0), P = params.EF_PIN;
    let eu = EF_BASE, ed = EF_BASE, el = EF_BASE, er = EF_BASE;
    const rootward = ((m.stele[i]||m.endo[i]) && !m.tip[i]) || (m.cortex[i] && r>=r_mer && !m.tip[i]);
    const shootward = (m.cortex[i] && r<r_mer) || (m.epi[i] && r<r_mer) || (m.epi[i] && m.tip[i]) || m.lrc[i];
    if (rootward) ed = P;
    if (shootward) eu = P;
    if (m.epi[i] && r >= r_mer) { if (c < c0) er = P; else if (c > c0) el = P; }
    if (m.col[i]) { eu = ed = el = er = P; if (params.gravity) { el = P*1.6; er = P*0.4; } }
    ef_up[i]=eu; ef_dn[i]=ed; ef_lf[i]=el; ef_rt[i]=er;
    pin_in[i] = (m.epi[i]||m.lrc[i]||m.col[i]||m.qc[i]||(m.cortex[i]&&m.tip[i])) ? params.IN_AUX1 : IN_BASE;
    prod[i] = m.qc[i] ? params.P_qc : (m.col[i] ? P_col : P_base);
  }
  for (let c = c0-1; c <= c0+1; c++) prod[c] += S_shoot;     // shoot supply at top of stele
}

let a = new Float64Array(N);
function reset() { a = new Float64Array(N); for (let i=0;i<N;i++) if (active[i]) a[i]=0.1; }

function step() {
  const da = new Float64Array(N);
  for (let i = 0; i < N; i++) if (active[i]) da[i] = prod[i] - kDecay * a[i];
  for (let r = 0; r < R; r++) for (let c = 0; c < C-1; c++) {
    const i=r*C+c, j=i+1; if (!(active[i]&&active[j])) continue;
    const W = (ef_rt[i]*a[i] + ef_lf[j]*a[j]) / (pin_in[i]+pin_in[j]);
    da[i] += g*(pin_in[i]*W - ef_rt[i]*a[i]);
    da[j] += g*(pin_in[j]*W - ef_lf[j]*a[j]);
  }
  for (let r = 0; r < R-1; r++) for (let c = 0; c < C; c++) {
    const i=r*C+c, j=i+C; if (!(active[i]&&active[j])) continue;
    const W = (ef_dn[i]*a[i] + ef_up[j]*a[j]) / (pin_in[i]+pin_in[j]);
    da[i] += g*(pin_in[i]*W - ef_dn[i]*a[i]);
    da[j] += g*(pin_in[j]*W - ef_up[j]*a[j]);
  }
  for (let r = 0; r < R; r++) for (let c = 0; c < C; c++) {     // boundary efflux
    const i=r*C+c; if (!active[i]) continue;
    if (c===0   || !active[i-1]) da[i] -= g*ef_lf[i]*a[i];
    if (c===C-1 || !active[i+1]) da[i] -= g*ef_rt[i]*a[i];
    if (r===0   || !active[i-C]) da[i] -= g*ef_up[i]*a[i];
    if (r===R-1 || !active[i+C]) da[i] -= g*ef_dn[i]*a[i];
  }
  for (let i = 0; i < N; i++) if (active[i]) a[i] += dt * da[i];
}

// ---- magma colour map with gamma (auxin spans orders of magnitude) ----
const STOPS = [[0,0,4],[59,15,112],[140,41,129],[222,73,104],[254,159,109],[252,253,191]];
const GAMMA = 0.45;
function magma(t){ t=Math.pow(Math.max(0,Math.min(1,t)),GAMMA); const x=t*(STOPS.length-1), i=Math.floor(x), f=x-i;
  const A=STOPS[i], B=STOPS[Math.min(i+1,STOPS.length-1)];
  return `rgb(${A[0]+(B[0]-A[0])*f|0},${A[1]+(B[1]-A[1])*f|0},${A[2]+(B[2]-A[2])*f|0})`; }

const cv = document.getElementById('root'); const ctx = cv.getContext('2d');
const PX = 11; cv.width = C*PX; cv.height = R*PX;
function render(){
  let max = 1e-6, mi = -1;
  for (let i=0;i<N;i++) if (active[i] && a[i]>max){ max=a[i]; mi=i; }
  ctx.fillStyle = '#0b0f0a'; ctx.fillRect(0,0,cv.width,cv.height);
  for (let r=0;r<R;r++) for (let c=0;c<C;c++){ const i=r*C+c; if(!active[i])continue;
    ctx.fillStyle = magma(a[i]/max); ctx.fillRect(c*PX, r*PX, PX, PX); }
  if (mi>=0){ const r=(mi/C|0), c=mi%C; ctx.strokeStyle='#7CFF6B'; ctx.lineWidth=2;
    ctx.beginPath(); ctx.arc(c*PX+PX/2, r*PX+PX/2, PX*0.7, 0, 7); ctx.stroke();
    document.getElementById('maxinfo').textContent =
      `auxin max ${max.toFixed(1)} µM at row ${r}` + ((r>=r_qc-1 && Math.abs(c-c0)<=2) ? '  ✓ at the QC' : ''); }
}

let running = true;
function frame(){ if (running) for (let s=0;s<300;s++) step(); render(); requestAnimationFrame(frame); }

function bind(id, key){ const el=document.getElementById(id), out=document.getElementById(id+'v');
  const set=()=>{ params[key]=parseFloat(el.value); if(out)out.textContent=params[key]; rebuild(); };
  el.addEventListener('input', set); set(); }
bind('pin','EF_PIN'); bind('aux1','IN_AUX1'); bind('prodqc','P_qc'); bind('decay','k');
document.getElementById('gravity').addEventListener('change', e=>{ params.gravity=e.target.checked; rebuild(); });
document.getElementById('reset').addEventListener('click', reset);
document.getElementById('pause').addEventListener('click', e=>{ running=!running; e.target.textContent=running?'Pause':'Play'; });

rebuild(); reset(); frame();
