// Virtual Root — in-browser auxin transport model (v2) + teaching/research UI:
// anatomy labels, reflux flux arrows, mutant/gravitropism presets, animated bending,
// PNG export, DII-VENUS compare view, hover-to-probe, µM/time/steady readouts, shareable links.

const R = 52, C = 17, c0 = 8, N = R * C;
const r_qc = 46, r_mer = 30, CAP = 4.7;
const L = 10.0, g = 1.0 / L, dt = 0.04;
const EF_BASE = 1, IN_BASE = 20;
const S_shoot = 0.30, P_col = 0.012, P_base = 0.005;

// ---- static tissue masks ----
const active = new Uint8Array(N);
const m = { stele:U(), endo:U(), cortex:U(), epi:U(), qc:U(), col:U(), lrc:U(), tip:U() };
function U(){ return new Uint8Array(N); }
for (let r = 0; r < R; r++) for (let c = 0; c < C; c++) {
  const i = r*C+c, rad = Math.abs(c - c0);
  active[i] = (r < r_qc) ? (rad <= 4 ? 1 : 0)
                         : (rad <= Math.sqrt(Math.max(0, CAP*CAP - (r-r_qc)*(r-r_qc))) ? 1 : 0);
  if (!active[i]) continue;
  m.stele[i]=rad<=1?1:0; m.endo[i]=rad===2?1:0; m.cortex[i]=rad===3?1:0; m.epi[i]=rad===4?1:0;
  m.tip[i]=r>=r_qc?1:0; m.qc[i]=(r===r_qc&&rad<=1)?1:0;
  m.col[i]=(r>r_qc&&rad<=2)?1:0; m.lrc[i]=(r>=r_qc-2&&rad>=3)?1:0;
}

// ---- parameter fields ----
const ef_up=new Float64Array(N), ef_dn=new Float64Array(N), ef_lf=new Float64Array(N), ef_rt=new Float64Array(N);
const pin_in=new Float64Array(N), prod=new Float64Array(N);
let kDecay = 0.0012;
const params = { EF_PIN:16, IN_AUX1:42, P_qc:0.025, k:0.0012, gravity:false, pin2ko:false, aux1ko:false,
                 labels:true, arrows:true, dii:false };

function rebuild(){
  kDecay = params.k; converged=false;
  for (let r=0;r<R;r++) for (let c=0;c<C;c++){
    const i=r*C+c;
    if(!active[i]){ ef_up[i]=ef_dn[i]=ef_lf[i]=ef_rt[i]=pin_in[i]=prod[i]=0; continue; }
    const rad=Math.abs(c-c0), P=params.EF_PIN;
    let eu=EF_BASE, ed=EF_BASE, el=EF_BASE, er=EF_BASE;
    const rootward=((m.stele[i]||m.endo[i])&&!m.tip[i])||(m.cortex[i]&&r>=r_mer&&!m.tip[i]);
    const shootward=(m.cortex[i]&&r<r_mer)||(m.epi[i]&&r<r_mer)||(m.epi[i]&&m.tip[i])||m.lrc[i];
    if(rootward) ed=P;
    if(shootward && !params.pin2ko) eu=P;
    if(m.epi[i]&&r>=r_mer){ if(c<c0) er=P; else if(c>c0) el=P; }
    if(m.col[i]){ eu=ed=el=er=P; if(params.gravity){ el=P*1.6; er=P*0.4; } }
    ef_up[i]=eu; ef_dn[i]=ed; ef_lf[i]=el; ef_rt[i]=er;
    const hasAux1=(m.epi[i]||m.lrc[i]||m.col[i]||m.qc[i]||(m.cortex[i]&&m.tip[i])) && !params.aux1ko;
    pin_in[i]=hasAux1?params.IN_AUX1:IN_BASE;
    prod[i]=m.qc[i]?params.P_qc:(m.col[i]?P_col:P_base);
  }
  for(let c=c0-1;c<=c0+1;c++) prod[c]+=S_shoot;
}

let a=new Float64Array(N), simTime=0, converged=false;
function reset(){ a=new Float64Array(N); for(let i=0;i<N;i++) if(active[i]) a[i]=0.1; simTime=0; converged=false; }

function step(){
  const da=new Float64Array(N);
  for(let i=0;i<N;i++) if(active[i]) da[i]=prod[i]-kDecay*a[i];
  for(let r=0;r<R;r++) for(let c=0;c<C-1;c++){ const i=r*C+c,j=i+1; if(!(active[i]&&active[j]))continue;
    const W=(ef_rt[i]*a[i]+ef_lf[j]*a[j])/(pin_in[i]+pin_in[j]);
    da[i]+=g*(pin_in[i]*W-ef_rt[i]*a[i]); da[j]+=g*(pin_in[j]*W-ef_lf[j]*a[j]); }
  for(let r=0;r<R-1;r++) for(let c=0;c<C;c++){ const i=r*C+c,j=i+C; if(!(active[i]&&active[j]))continue;
    const W=(ef_dn[i]*a[i]+ef_up[j]*a[j])/(pin_in[i]+pin_in[j]);
    da[i]+=g*(pin_in[i]*W-ef_dn[i]*a[i]); da[j]+=g*(pin_in[j]*W-ef_up[j]*a[j]); }
  for(let r=0;r<R;r++) for(let c=0;c<C;c++){ const i=r*C+c; if(!active[i])continue;
    if(c===0||!active[i-1]) da[i]-=g*ef_lf[i]*a[i];
    if(c===C-1||!active[i+1]) da[i]-=g*ef_rt[i]*a[i];
    if(r===0||!active[i-C]) da[i]-=g*ef_up[i]*a[i];
    if(r===R-1||!active[i+C]) da[i]-=g*ef_dn[i]*a[i]; }
  let md=0; for(let i=0;i<N;i++) if(active[i]){ const d=dt*da[i]; a[i]+=d; if(Math.abs(d)>md)md=Math.abs(d); }
  return md;
}

// ---- colour maps (gamma; auxin spans orders of magnitude) ----
const MAGMA=[[0,0,4],[59,15,112],[140,41,129],[222,73,104],[254,159,109],[252,253,191]];
const VENUS=[[4,6,2],[18,55,20],[70,140,45],[170,200,60],[250,252,205]]; // YFP-like, for DII-VENUS view
const GAMMA=0.45;
function ramp(stops,t){ t=Math.pow(Math.max(0,Math.min(1,t)),GAMMA); const x=t*(stops.length-1),i=Math.floor(x),f=x-i;
  const A=stops[i],B=stops[Math.min(i+1,stops.length-1)];
  return `rgb(${A[0]+(B[0]-A[0])*f|0},${A[1]+(B[1]-A[1])*f|0},${A[2]+(B[2]-A[2])*f|0})`; }

// ---- canvas ----
const cv=document.getElementById('root'), ctx=cv.getContext('2d');
const PX=11, OX=64, OY=10, PADR=92, PADB=48;
cv.width=OX+C*PX+PADR; cv.height=OY+R*PX+PADB;
const cx=c=>OX+c*PX+PX/2, cy=r=>OY+r*PX+PX/2;

function arrow(x,y,dx,dy,len,col,lw){ const mg=Math.hypot(dx,dy)||1; dx/=mg; dy/=mg;
  const ex=x+dx*len, ey=y+dy*len, ang=Math.atan2(dy,dx), h=4.5;
  ctx.strokeStyle=col; ctx.fillStyle=col; ctx.lineWidth=lw;
  ctx.beginPath(); ctx.moveTo(x-dx*len,y-dy*len); ctx.lineTo(ex,ey); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(ex,ey);
  ctx.lineTo(ex-h*Math.cos(ang-0.5),ey-h*Math.sin(ang-0.5));
  ctx.lineTo(ex-h*Math.cos(ang+0.5),ey-h*Math.sin(ang+0.5));
  ctx.closePath(); ctx.fill(); }

function label(txt,tx,ty,tcol,gx,gy){ ctx.strokeStyle='#889'; ctx.lineWidth=0.7;
  ctx.beginPath(); ctx.moveTo(tx<OX?tx+2:tx-2,ty); ctx.lineTo(gx,gy); ctx.stroke();
  ctx.fillStyle=tcol; ctx.font='11px system-ui,sans-serif';
  ctx.textAlign=tx<OX?'right':'left'; ctx.textBaseline='middle'; ctx.fillText(txt,tx,ty); }

let curMax=1;
function render(){
  let max=1e-6, mi=-1;
  for(let i=0;i<N;i++) if(active[i]&&a[i]>max){ max=a[i]; mi=i; }
  curMax=max;
  ctx.fillStyle='#0b0f0a'; ctx.fillRect(0,0,cv.width,cv.height);
  const stops = params.dii ? VENUS : MAGMA;
  for(let r=0;r<R;r++) for(let c=0;c<C;c++){ const i=r*C+c; if(!active[i])continue;
    const t=a[i]/max; ctx.fillStyle=ramp(stops, params.dii ? 1-t : t); ctx.fillRect(OX+c*PX,OY+r*PX,PX,PX); }

  if(params.arrows){
    for(let r=2;r<R;r+=3) for(let c=1;c<C;c+=2){ const i=r*C+c; if(!active[i])continue;
      const vx=(ef_rt[i]-ef_lf[i])*a[i], vy=(ef_dn[i]-ef_up[i])*a[i];
      if(Math.abs(vx)+Math.abs(vy)>0.4) arrow(cx(c),cy(r),vx,vy,4.2,'rgba(255,255,255,0.5)',1.1); }
  }
  if(mi>=0){ const r=mi/C|0,c=mi%C; ctx.strokeStyle='#7CFF6B'; ctx.lineWidth=2;
    ctx.beginPath(); ctx.arc(cx(c),cy(r),PX*0.7,0,7); ctx.stroke(); }
  if(params.labels){
    label('Elongation zone',OX-10,cy(14),'#dfe','#fff',cx(4)-PX/2,cy(14));
    label('Meristem',OX-10,cy(32),'#dfe',OX+(c0-4)*PX,cy(32));
    label('QC',OX-10,cy(r_qc),'#7CFF6B',cx(c0)-PX,cy(r_qc));
    label('Columella / cap',OX-10,cy(r_qc+2),'#dfe',OX+2*PX,cy(r_qc+2));
    const RX=OX+C*PX+10;
    label('Stele',RX,cy(20),'#dfe',cx(c0),cy(20));
    label('Cortex',RX,cy(24),'#dfe',cx(c0+3),cy(24));
    label('Epidermis',RX,cy(28),'#dfe',cx(c0+4),cy(28));
  }
  let lo=0,hi=0;
  for(let r=r_mer;r<r_qc-2;r++) for(let c=0;c<C;c++){ const i=r*C+c; if(!active[i])continue;
    if(c<c0)lo+=a[i]; else if(c>c0)hi+=a[i]; }
  const asym=(lo-hi)/(lo+hi+1e-9), by=OY+R*PX+22;
  if(Math.abs(asym)>0.03){ const dir=asym>0?-1:1;
    arrow(cx(c0),by,dir,0,10+40*Math.min(1,Math.abs(asym)*4),'#39c6ff',2.4);
    ctx.fillStyle='#39c6ff'; ctx.font='11px system-ui'; ctx.textAlign='center'; ctx.textBaseline='top';
    ctx.fillText('root bends toward the lower side',cx(c0),by+8); }

  const atQC = mi>=0 && (mi/C|0)>=r_qc-1 && Math.abs(mi%C-c0)<=2;
  document.getElementById('maxinfo').textContent =
    (mi>=0?`auxin max ${max.toFixed(1)} µM at the ${atQC?'QC ✓':'row '+(mi/C|0)}`:'') +
    (params.gravity?` · gravistimulated (${(Math.abs(asym)*100).toFixed(0)}%)`:'');
  document.getElementById('legmax').textContent = max.toFixed(1)+' µM';
  document.getElementById('clock').textContent =
    `t = ${(simTime/60).toFixed(1)} min` + (converged?' · steady state ✓':' · settling…') +
    (params.dii?' · DII-VENUS view (bright = low auxin)':'');
}

// ---- animation ----
let running=true;
function frame(){ if(running && !converged){ let md=0; for(let s=0;s<300;s++){ const d=step(); if(d>md)md=d; simTime+=dt; } if(md<1e-6) converged=true; } render(); requestAnimationFrame(frame); }

// ---- probe (hover / touch to inspect a cell) ----
const probeEl=document.getElementById('probe');
function tissue(i){ return m.qc[i]?'QC':m.col[i]?'columella':m.lrc[i]?'root cap':m.epi[i]?'epidermis':m.cortex[i]?'cortex':m.endo[i]?'endodermis':'stele'; }
cv.addEventListener('pointermove', e=>{ const rect=cv.getBoundingClientRect();
  const c=Math.floor(((e.clientX-rect.left)*(cv.width/rect.width)-OX)/PX);
  const r=Math.floor(((e.clientY-rect.top)*(cv.height/rect.height)-OY)/PX);
  const i=r*C+c;
  if(r>=0&&r<R&&c>=0&&c<C&&active[i]){ const pins=[];
    if(ef_up[i]>EF_BASE)pins.push('↑'); if(ef_dn[i]>EF_BASE)pins.push('↓');
    if(ef_lf[i]>EF_BASE)pins.push('←'); if(ef_rt[i]>EF_BASE)pins.push('→');
    probeEl.textContent=`${tissue(i)} · auxin ${a[i].toFixed(2)} µM · PIN ${pins.join('')||'—'} · AUX1 ${pin_in[i]>IN_BASE?'yes':'no'}`;
  } else probeEl.textContent='hover / tap a cell to inspect it'; });

// ---- presets ----
const PRESETS={
  wt:{EF_PIN:16,IN_AUX1:42,P_qc:0.025,k:0.0012,gravity:false,pin2ko:false,aux1ko:false},
  grav:{EF_PIN:16,IN_AUX1:42,P_qc:0.025,k:0.0012,gravity:true,pin2ko:false,aux1ko:false},
  pin2:{EF_PIN:16,IN_AUX1:42,P_qc:0.025,k:0.0012,gravity:false,pin2ko:true,aux1ko:false},
  aux1:{EF_PIN:16,IN_AUX1:42,P_qc:0.025,k:0.0012,gravity:false,pin2ko:false,aux1ko:true},
  nodecay:{EF_PIN:16,IN_AUX1:42,P_qc:0.025,k:0.0005,gravity:false,pin2ko:false,aux1ko:false},
};
function syncDOM(){
  document.getElementById('pin').value=params.EF_PIN;   document.getElementById('pinv').textContent=params.EF_PIN;
  document.getElementById('aux1').value=params.IN_AUX1; document.getElementById('aux1v').textContent=params.IN_AUX1;
  document.getElementById('prodqc').value=params.P_qc;  document.getElementById('prodqcv').textContent=params.P_qc;
  document.getElementById('decay').value=params.k;      document.getElementById('decayv').textContent=params.k;
  document.getElementById('gravity').checked=params.gravity; document.getElementById('dii').checked=params.dii;
}
function setPreset(name){ Object.assign(params,PRESETS[name]); syncDOM(); rebuild(); reset(); }

// ---- shareable URL state ----
function encodeState(){ const p=params;
  return `p=${p.EF_PIN}&a=${p.IN_AUX1}&q=${p.P_qc}&k=${p.k}&g=${+p.gravity}&m2=${+p.pin2ko}&m1=${+p.aux1ko}&d=${+p.dii}`; }
function applyState(str){ const u=new URLSearchParams(str); if(!u.has('p')) return false;
  const n=(k,d)=>u.has(k)?+u.get(k):d;
  params.EF_PIN=n('p',16); params.IN_AUX1=n('a',42); params.P_qc=n('q',0.025); params.k=n('k',0.0012);
  params.gravity=n('g',0)===1; params.pin2ko=n('m2',0)===1; params.aux1ko=n('m1',0)===1; params.dii=n('d',0)===1;
  syncDOM(); return true; }

// ---- UI ----
function bind(id,key){ const el=document.getElementById(id),out=document.getElementById(id+'v');
  const set=()=>{ params[key]=parseFloat(el.value); if(out)out.textContent=params[key]; rebuild(); };
  el.addEventListener('input',set); set(); }
bind('pin','EF_PIN'); bind('aux1','IN_AUX1'); bind('prodqc','P_qc'); bind('decay','k');
document.getElementById('gravity').addEventListener('change',e=>{params.gravity=e.target.checked;rebuild();});
document.getElementById('dii').addEventListener('change',e=>{params.dii=e.target.checked;});
document.getElementById('labels').addEventListener('change',e=>params.labels=e.target.checked);
document.getElementById('arrows').addEventListener('change',e=>params.arrows=e.target.checked);
document.getElementById('reset').addEventListener('click',reset);
document.getElementById('pause').addEventListener('click',e=>{running=!running;e.target.textContent=running?'Pause':'Play';});
document.getElementById('png').addEventListener('click',()=>{ const l=document.createElement('a');
  l.download='virtual-root.png'; l.href=cv.toDataURL('image/png'); l.click(); });
document.getElementById('share').addEventListener('click',e=>{ location.hash=encodeState();
  navigator.clipboard && navigator.clipboard.writeText(location.href);
  e.target.textContent='Link copied ✓'; setTimeout(()=>e.target.textContent='Copy share link',1500); });
document.querySelectorAll('[data-preset]').forEach(b=>b.addEventListener('click',()=>setPreset(b.dataset.preset)));

applyState(location.hash.slice(1)) || applyState(location.search.slice(1));
rebuild(); reset(); frame();
