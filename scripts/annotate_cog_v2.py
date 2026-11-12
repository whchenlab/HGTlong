"""Conservative COG homology annotation, not eggNOG orthology inference."""
import argparse,csv,json,re
from pathlib import Path
from collections import defaultdict,Counter
p=argparse.ArgumentParser();p.add_argument('--genes',required=True);p.add_argument('--hits',required=True);p.add_argument('--db',required=True);p.add_argument('--out',required=True);args=p.parse_args()
db=Path(args.db);defs={}
for line in open(db/'cog-24.def.tab'):
 f=line.rstrip('\n').split('\t');defs[f[0]]=dict(category=f[1],description=f[2],name=f[3] if len(f)>3 else '')
hits=defaultdict(list);needed=set()
for line in open(args.hits):
 f=line.rstrip().split('\t');q,s=f[:2];fields=dict(subject=s,identity=float(f[2]),qstart=int(f[4]),qend=int(f[5]),sstart=int(f[6]),send=int(f[7]),evalue=float(f[8]),bitscore=float(f[9]),qlen=int(f[10]),slen=int(f[11]))
 if fields['identity']<30 or fields['evalue']>1e-10 or (fields['qend']-fields['qstart']+1)/fields['qlen']<.7 or (fields['send']-fields['sstart']+1)/fields['slen']<.7:continue
 hits[q].append(fields);needed.add(s)
footprints=defaultdict(list)
for f in csv.reader(open(db/'cog-24.cog.csv')):
 if len(f)>=9 and f[2] in needed:
  spans=[tuple(map(int,x.split('-'))) for x in f[4].split('=')]
  entry=(f[6],tuple(spans),f[8])
  if entry not in footprints[f[2]]:footprints[f[2]].append(entry)
rows=[]
for row in csv.DictReader(open(args.genes),delimiter='\t'):
 h=hits.get(row['gene_id'],[]);best=max((x['bitscore'] for x in h),default=0);top=[x for x in h if x['bitscore']>=.95*best];sets=[]
 for x in top:
  cogs=set()
  for cog,spans,memberclass in footprints.get(x['subject'],[]):
   covered=sum(max(0,min(b,x['send'])-max(a,x['sstart'])+1) for a,b in spans)
   total=sum(b-a+1 for a,b in spans)
   if covered/total>=.7:cogs.add(cog)
  if cogs:sets.append(cogs)
 common=set.intersection(*sets) if sets else set()
 status='COG_homology_consensus' if common else 'conflicting_near_best_COGs' if sets else 'no_qualified_COG_hit'
 cats=''.join(sorted({c for cog in common for c in defs[cog]['category']}));descs='; '.join(defs[c]['description'] for c in sorted(common));names=';'.join(sorted({defs[c]['name'] for c in common if defs[c]['name']}))
 mobile=[];low=descs.lower()
 if re.search(r'transposase|transposition protein|insertion sequence',low):mobile.append('transposition')
 if re.search(r'integrase|integration host.*phage|site.specific recombinase|excisionase',low) and not re.search(r'\bxer[cd]\b|chromosom.*dimer|dif site',low+' '+names.lower()):mobile.append('integration_recombination')
 if re.search(r'conjugat|relaxase|mobilization protein|mobilisation protein',low):mobile.append('conjugation_mobilization')
 if re.search(r'phage.*(capsid|tail|portal|terminase)|(?:capsid|terminase).*phage|prophage',low):mobile.append('phage_associated')
 if 'X' in cats and not mobile:mobile.append('mobilome_COG_X')
 strongest=max(top,key=lambda x:x['bitscore']) if top else {}
 row.update(cog_ids=';'.join(sorted(common)),cog_categories=cats,gene_name=names,description=descs,annotation_status=status,mobile_types=';'.join(sorted(set(mobile))),best_subject=strongest.get('subject',''),best_identity=strongest.get('identity',''),best_evalue=strongest.get('evalue',''),best_bitscore=best or '',annotation_source='NCBI_COG2024_DIAMOND_homology')
 rows.append(row)
with open(args.out,'w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
Path(args.out+'.summary.json').write_text(json.dumps(dict(genes=len(rows),status=dict(Counter(r['annotation_status'] for r in rows)),mobile_genes=sum(bool(r['mobile_types']) for r in rows),parameters=dict(identity=30,evalue=1e-10,query_coverage=.7,subject_coverage=.7,COG_footprint_coverage=.7,near_best_bitscore_fraction=.95),interpretation='Provisional protein homology annotations; mobile-associated genes do not delineate complete mobile elements.'),indent=2)+'\n')
print('Annotated',len(rows),'genes',flush=True)
