"""Standalone, offline read-evidence browser. All text is escaped by DOM APIs."""
import json
from datetime import datetime,timezone,timedelta
from pathlib import Path
from collections import defaultdict


def write_report(events,members,path,context=None):
    tier={'read_supported':0,'partial':1,'ambiguous':2}
    events=sorted(events,key=lambda e:(-int(e.get('screened_candidate',0)),tier.get(e['status'],3),
        e.get('qc_status')=='repeat_associated',
        -min(int(e['n_'+k]) for k in ('a0','a1','b0','b1')),-int(e['molecules']),e['event_id']))
    evidence=defaultdict(list)
    by_side=defaultdict(lambda:defaultdict(list))
    seen=defaultdict(set)
    fields=('read','resident','resident_start','resident_end','resident_identity',
            'fragment_identity','ambiguous_fragment','assignment_ambiguous','a0','a1','b0','b1')
    for row in sorted(members,key=lambda r:(int(r.get('assignment_ambiguous',0)),int(r['ambiguous_fragment']))):
        key=(row['read'],row['resident'])
        eid=row['event_id']
        side=row['resident']
        if key not in seen[eid] and len(by_side[eid][side])<100:
            by_side[eid][side].append({k:row.get(k,0) for k in fields})
        seen[eid].add(key)
    for eid,sides in by_side.items():
        first=sides['A'][:50]+sides['B'][:50]
        remaining=sides['A'][50:]+sides['B'][50:]
        evidence[eid]=(first+remaining)[:100]
    data=json.dumps(dict(context=context,events=events,evidence=evidence,evidence_counts={k:len(v) for k,v in seen.items()},
                        generated_at=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')),
                    separators=(',',':')).replace('<','\\u003c')
    template='''<!doctype html><html><head><meta charset="utf-8"><title>HGTlong · Read evidence</title>
<style>body{font:15px system-ui;background:#f4f6fa;color:#17283a;margin:0}header,main{padding:24px 32px}header{background:#142f43;color:white}h1{margin:0}h2{font-size:20px}.cards{display:flex;gap:16px;margin:18px 0}.card,section{background:white;padding:18px;border-radius:10px;box-shadow:0 1px 4px #ccd5de}.card b{display:block;font-size:28px;color:#087d88}input{padding:10px;width:400px;max-width:80%}table{border-collapse:collapse;width:100%;font-size:12px}td,th{padding:8px;border-bottom:1px solid #dce3e9;text-align:left}th{background:#edf1f5;position:sticky;top:0}tr[data-id]{cursor:pointer}tr[data-id]:hover{background:#dff3f4}.scroll{overflow:auto;max-height:420px}.note{color:#57697a}svg{width:100%;min-width:800px}button{padding:8px;cursor:pointer}.tag{font-family:monospace}#hover-tip{position:fixed;z-index:9999;max-width:540px;padding:12px 14px;border-radius:8px;background:#142f43;color:white;box-shadow:0 4px 18px #0005;white-space:pre-line;font:13px/1.5 system-ui;pointer-events:none}#hover-tip[hidden]{display:none}[data-tooltip]{cursor:help}.read-alignment.read-linked{stroke:#b51c91;stroke-width:7;opacity:1}.hdmi-final{color:#986100;font-weight:700}</style></head>
<body><header><h1>HGTlong</h1><p>ONT read-resolved candidate loci · Research release candidate 0.1.3rc3</p><small id="stamp"></small></header><main>
<p>Candidate loci are defined from individual ONT molecules. Read support does not establish evolutionary transfer or its direction. HDMI is used only in a separate concordance analysis.</p>
<div class="cards" id="counts"></div><section><h2>Candidate loci</h2><input id="filter" placeholder="Search event, MAG, contig, or status"><p class="note">Select a row to inspect its molecules. All coordinates are 0-based, half-open.</p><p id="matchcount" class="note"></p><div class="scroll"><table><thead id="head"></thead><tbody id="rows"></tbody></table></div></section>
<section style="margin-top:20px" id="detail"><h2>Select a candidate</h2></section>
<p class="note">Repeat screening is shown per candidate; supplied masks may not cover all repeat classes. Unique read names approximate molecules; split/duplex metadata is not available. A complete fragment is required to seed resolved endpoints. All placements are retained in event_evidence.tsv.gz. See README.md for calling semantics and limitations.</p>
</main><script id="data" type="application/json">DATA</script><script>
const data=JSON.parse(document.getElementById('data').textContent), cols=['event_id','status','magA','magB','molecules','n_a0','n_a1','n_b0','n_b1','length_A','length_B','length_ratio'];
for(const e of data.events){e.length_A=Number(e.a1)-Number(e.a0);e.length_B=Number(e.b1)-Number(e.b0);}
if(data.events.some(e=>'hdmi_category' in e))cols.push('hdmi_category');
if(data.events.some(e=>'screened_candidate' in e))cols.splice(2,0,'screened_candidate');
document.getElementById('stamp').textContent='Generated '+data.generated_at+' (UTC+08:00)';
function el(tag,text,parent){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(parent)parent.append(e);return e;}
const counts={};if(data.events.some(e=>'screened_candidate' in e))counts.screened_candidates=data.events.filter(e=>Number(e.screened_candidate)===1).length;for(const e of data.events)counts[e.status]=(counts[e.status]||0)+1;for(const [k,v] of Object.entries(counts)){const c=el('div',undefined,document.getElementById('counts'));c.className='card';el('b',v,c);el('span',k,c);}
const h=el('tr',undefined,document.getElementById('head'));for(const c of cols)el('th',c,h);
function render(){const body=document.getElementById('rows');body.replaceChildren();const q=document.getElementById('filter').value.toLowerCase();const matches=data.events.filter(e=>JSON.stringify(e).toLowerCase().includes(q));document.getElementById('matchcount').textContent='Showing '+Math.min(matches.length,500)+' of '+matches.length+' matches'+(matches.length>500?' — narrow the search to inspect others.':'.');for(const e of matches.slice(0,500)){const r=el('tr',undefined,body);r.dataset.id=e.event_id;for(const c of cols){const full=String(e[c]??'not assessed'),text=c==='magA'||c==='magB'?(full.match(/cleanbin_\\d+$/)||[full])[0]:full;el('td',text,r).title=full;}r.onclick=()=>detail(e);}}
const tip=el('div',undefined,document.body);tip.id='hover-tip';tip.setAttribute('role','tooltip');tip.hidden=true;
function positionTip(x,y){tip.style.left=Math.max(8,Math.min(x+14,innerWidth-tip.offsetWidth-8))+'px';tip.style.top=Math.max(8,Math.min(y+14,innerHeight-tip.offsetHeight-8))+'px';}
function linkRead(target){const id=target?.dataset.readId;for(const line of document.querySelectorAll(".read-alignment"))line.classList.toggle("read-linked",id!==undefined&&line.dataset.readId===id);}
function hideTip(){tip.hidden=true;linkRead(null);}
function showTip(target,x,y){linkRead(target);tip.textContent=target.dataset.tooltip;tip.hidden=false;positionTip(x,y);}
document.addEventListener('pointerover',event=>{const target=event.target.closest?.('[data-tooltip]');if(target)showTip(target,event.clientX,event.clientY);});
document.addEventListener('pointermove',event=>{if(!tip.hidden)positionTip(event.clientX,event.clientY);});
document.addEventListener('pointerout',event=>{if(event.target.closest?.('[data-tooltip]'))hideTip();});
document.addEventListener('focusin',event=>{const target=event.target.closest?.('[data-tooltip]');if(target){const r=target.getBoundingClientRect();showTip(target,r.left,r.bottom);}});
document.addEventListener('focusout',hideTip);document.addEventListener('keydown',event=>{if(event.key==='Escape')hideTip();});
function svgEl(tag,attrs,parent,text){if(tag==='title'){parent.dataset.tooltip=text;parent.setAttribute('tabindex','0');parent.setAttribute('aria-label',text);return parent;}const x=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))x.setAttribute(k,v);if(text!==undefined)x.textContent=text;parent.append(x);if(tag==='title'){parent.dataset.tooltip=text;parent.setAttribute('tabindex','0');parent.setAttribute('aria-label',text);}return x;}
let flank=5000,hdmiScope=data.context?.scope_counts?.final_retained!==undefined?"final":"coordinate";
function contextTracks(e,box){
 const ctx=data.context;if(!ctx)return;
 const controls=el('p','Genomic neighborhood: ',box),choice=el('select',undefined,controls);
 for(const n of [1000,5000,10000]){const o=el('option',(n/1000)+' kb flank',choice);o.value=n;o.selected=n===flank;}
 choice.onchange=()=>{flank=Number(choice.value);detail(e);};
 el('span',' · HDMI comparator: ',controls);const scope=el('select',undefined,controls);for(const [value,label] of [...(ctx.scope_counts?.final_retained!==undefined?[['final','Final retained ('+ctx.scope_counts.final_retained+')']]:[]),...(ctx.scope_counts?.software_retained!==undefined?[['software','HDMI software summary ('+ctx.scope_counts.software_retained+')']]:[]),['coordinate','Coordinate-pass ('+(ctx.scope_counts?.coordinate_pass??'?')+')'],['raw','Raw ('+(ctx.scope_counts?.raw??'?')+')'],['reviewed','Reviewed subset ('+(ctx.scope_counts?.reviewed??'?')+')']]){const o=el('option',label,scope);o.value=value;o.selected=value===hdmiScope;}scope.onchange=()=>{hdmiScope=scope.value;detail(e);};
 el('p','Arrows: predicted CDS in native reference orientation; color: COG category; red outline: mobile-associated annotation. Gold: final retained HDMI candidates; purple: other HDMI intervals. Tracks can include events with other partners. These are nearby annotations, not demonstrated transfer mechanisms.',box);
 const sides=[['A',20,Number(e.a0),Number(e.a1),e.ctgA],['B',580,Number(e.b0),Number(e.b1),e.ctgB]],prepared=[];
 for(const [side,left,start,end,ctg] of sides){
  const lo=Math.max(0,start-flank),hi=Math.min(ctx.lengths[ctg]||end+flank,end+flank);
  const genes=(ctx.genes[ctg]||[]).filter(g=>g.end>lo&&g.start<hi),hdmi=(ctx.hdmi[ctg]||[]).filter(h=>h.end>lo&&h.start<hi&&(hdmiScope==='raw'||(hdmiScope==='final'?h.final_retained:hdmiScope==='software'?h.software_retained:hdmiScope==='reviewed'?h.reviewed:h.coordinate_pass)));
  const laneEnds=[];for(const g of genes){let lane=laneEnds.findIndex(v=>v<=g.start);if(lane<0)lane=laneEnds.length;laneEnds[lane]=g.end;g.lane=lane;}
  prepared.push({side,left,start,end,ctg,lo,hi,genes,hdmi,lanes:laneEnds.length});
 }
 const height=Math.max(...prepared.map(s=>110+s.lanes*26+s.hdmi.length*18)),wrap=el('div',undefined,box);wrap.style.overflow='auto';const svg=svgEl('svg',{viewBox:'0 0 1120 '+height},wrap);
 for(const s of prepared){const x=p=>s.left+520*Math.max(0,Math.min(1,(p-s.lo)/(s.hi-s.lo)));
  svgEl('text',{x:s.left,y:18,'font-size':12},svg,s.side+' · '+s.lo.toLocaleString()+'–'+s.hi.toLocaleString());
  svgEl('rect',{x:x(s.start),y:26,width:Math.max(1,x(s.end)-x(s.start)),height:height-40,fill:'#e5eef4'},svg);
  for(const p of [s.start,s.end])svgEl('line',{x1:x(p),x2:x(p),y1:26,y2:height-12,stroke:'#7a8b9b','stroke-dasharray':'4 3'},svg);
  for(const g of s.genes){const a=x(g.start),b=x(g.end),y=42+g.lane*26,tip=Math.min(7,(b-a)/2),points=g.strand==='+'?[[a,y-6],[b-tip,y-6],[b,y],[b-tip,y+6],[a,y+6]]:[[b,y-6],[a+tip,y-6],[a,y],[a+tip,y+6],[b,y+6]];
   const category=g.cog_categories[0],color=(ctx.categories[category]||{}).color||'#bbc4cc';
   const shape=svgEl('polygon',{points:points.map(p=>p.join(',')).join(' '),fill:color,stroke:g.mobile_types?'#c73535':'#586a7c','stroke-width':g.mobile_types?2:0.6},svg);
   const label=g.gene_name||g.cog_ids||'unassigned';if(b-a>42)svgEl('text',{x:(a+b)/2,y:y-9,'font-size':8,'text-anchor':'middle'},svg,label.slice(0,18));
   svgEl('title',{},shape,g.gene_id+'\\n'+g.start+'–'+g.end+' ('+g.strand+')\\nCOG: '+(g.cog_ids||'unassigned')+' | categories: '+(g.cog_categories||'unassigned')+'\\n'+(g.description||'No qualified COG annotation')+'\\nMobile-associated: '+(g.mobile_types||'not identified in these annotations')+'\\n'+g.annotation_status);
  }
  const base=65+s.lanes*26;svgEl('text',{x:s.left,y:base,'font-size':10},svg,'Nearby HDMI ('+s.hdmi.length+')');
  s.hdmi.forEach((h,i)=>{const y=base+14+i*18,a=x(h.start),b=x(h.end),line=svgEl('line',{x1:a,x2:b,y1:y,y2:y,stroke:h.final_retained?'#c18400':'#7852a9','stroke-width':h.final_retained?7:5},svg);svgEl('title',{},line,h.id+(h.final_retained?' · FINAL RETAINED':' · not in final set')+(h.partner===(s.side==='A'?e.ctgB:e.ctgA)?' · same contig pair':' · other partner')+' · '+h.start+'–'+h.end+'\\nPartner: '+h.partner+'\\nOrientation: '+h.orientation+'\\nRetention: '+(h.retention_reason||'not supplied'));svgEl('text',{x:a,y:y-4,'font-size':8},svg,(h.final_retained?'★ ':'')+h.id);});
 }
 const legend=el('details',undefined,box);el('summary','COG color legend',legend);for(const c of [...new Set(prepared.flatMap(s=>s.genes.map(g=>g.cog_categories[0])).filter(Boolean))].sort()){const label=el('p',c+' · '+(ctx.categories[c]?.name||'unassigned'),legend);label.style.borderLeft='8px solid '+(ctx.categories[c]?.color||'#bbc4cc');label.style.paddingLeft='8px';}
 const details=el('details',undefined,box);el('summary','Genes overlapping the block and mobile-associated genes in the displayed flanks',details);const table=el('table',undefined,details),head=el('tr',undefined,table);for(const k of ['Side','Gene / coordinates','COG categories','Function','Mobile association'])el('th',k,head);
 for(const s of prepared)for(const g of s.genes){const inside=g.end>s.start&&g.start<s.end;if(!inside&&!g.mobile_types)continue;const r=el('tr',undefined,table);for(const t of [s.side,g.gene_id+' · '+g.start+'–'+g.end+(inside?' · block overlap':' · flank'),g.cog_categories||'unassigned',g.description||'No qualified COG annotation',g.mobile_types||'none identified'])el('td',t,r);}
 el('p',ctx.note,box).className='note';
}

function detail(e){tip.hidden=true;const box=document.getElementById('detail');box.replaceChildren();el('h2',e.event_id+' · '+e.status,box);el('p',e.ctgA+':'+e.a0+'–'+e.a1+' ↔ '+e.ctgB+':'+e.b0+'–'+e.b1+' | orientation '+e.orientation+' | block lengths A/B: '+e.length_A+' / '+e.length_B+' bp',box);el('p','Four-boundary support: '+[e.n_a0,e.n_a1,e.n_b0,e.n_b1].join(' / ')+' | ambiguous molecules: '+e.ambiguous_molecules+' | ambiguously assigned molecules: '+(e.unassigned_molecules||0)+' | repeat screen ('+(e.mask_classes||'unspecified masks')+'): '+e.repeat_screen,box);if(e.sequence_qc)el('p','Block sequence QC: '+e.sequence_qc+' | aligned identity: '+e.block_identity+' | aligned fractions A/B: '+e.aligned_fraction_A+' / '+e.aligned_fraction_B,box);
if('screened_candidate' in e)el('p','Screened candidate: '+(Number(e.screened_candidate)===1?'yes':'no')+' | Failed criteria: '+(e.screen_fail_reasons||'none')+'. Screening requires ≥3 reads at each boundary, ≥99% aligned block identity, ≥95% coverage on each block, and supplied annotation QC pass. rRNA-only masks do not assess other repeat classes.',box);
if(e.hdmi_category){el('p','HDMI comparison: '+e.hdmi_category,box).style.fontWeight='700';el('p',e.hdmi_reason,box);}
contextTracks(e,box);
const shown=data.evidence[e.event_id]||[];el('p','Showing '+shown.length+' of '+(data.evidence_counts[e.event_id]||0)+' molecule placements, balanced between resident sides where possible. Teal: resident alignment; orange: partner fragment. Hover or focus a read to highlight the same molecule in both blocks (magenta) and show its identities.',box);
const scroll=el('div',undefined,box);scroll.style.overflow='auto';const svg=svgEl('svg',{viewBox:'0 0 1120 '+(95+shown.length*14)},scroll);
for(const [side,left,start,end,ctg] of [['A',20,e.a0,e.a1,e.ctgA],['B',580,e.b0,e.b1,e.ctgB]]){const lo=Math.max(0,start-flank),hi=Math.min(data.context?.lengths[ctg]||end+flank,end+flank),x=p=>left+Math.max(0,Math.min(1,(p-lo)/(hi-lo)))*520;svgEl('text',{x:left,y:16,'font-size':11},svg,side+' · '+ctg);svgEl('text',{x:left,y:36,'font-size':14,'font-weight':'700',fill:'#087d88','class':'read-block-length'},svg,'HGT block '+side+': '+(end-start).toLocaleString()+' bp');svgEl('rect',{x:x(start),y:45,width:x(end)-x(start),height:35+shown.length*14,fill:'#e5eef4'},svg);for(const p of [start,end]){svgEl('line',{x1:x(p),x2:x(p),y1:45,y2:80+shown.length*14,stroke:'#889aaa','stroke-dasharray':'3 3'},svg);svgEl('text',{x:x(p)+(p===start?-2:2),y:58,'font-size':10,'text-anchor':p===start?'end':'start'},svg,p);}
shown.forEach((r,i)=>{const resident=r.resident===side,s=resident?r.resident_start:r[side==='A'?'a0':'b0'],t=resident?r.resident_end:r[side==='A'?'a1':'b1'];const line=svgEl('line',{x1:x(s),x2:x(t),y1:75+i*14,y2:75+i*14,'class':'read-alignment','data-read-id':r.read,'data-block':side,stroke:resident?'#087d88':'#c77725','stroke-width':4,opacity:(Number(r.ambiguous_fragment)||Number(r.assignment_ambiguous)) ? 0.4 : 1},svg);svgEl('title',{},line,r.read+'\\nIdentity to A: '+(100*Number(r.resident==='A'?r.resident_identity:r.fragment_identity)).toFixed(3)+'%\\nIdentity to B: '+(100*Number(r.resident==='B'?r.resident_identity:r.fragment_identity)).toFixed(3)+'%\\nAlignment: '+(resident?'resident alignment':'partner fragment')+' '+s+'–'+t+'\\nNote: identities refer to each retained alignment, not the identical query interval');});}
el('p','Read-supported is a geometric evidence tier, not an HGT truth label. Transfer direction remains unresolved.',box);}
document.getElementById('filter').oninput=render;render();if(data.events.length)detail(data.events.find(e=>e.event_id===location.hash.slice(1))||data.events[0]);else document.getElementById('detail').textContent='No candidates met the current calling thresholds. Consult extraction statistics to distinguish low support from missing usable mappings.';
window.addEventListener('hashchange',()=>{const e=data.events.find(e=>e.event_id===location.hash.slice(1));if(e)detail(e);});
</script></body></html>'''
    Path(path).write_text(template.replace('DATA</script>',data+'</script>'))
