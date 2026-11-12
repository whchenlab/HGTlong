import csv,json,sys
from pathlib import Path

root=Path(sys.argv[1]); results=[]
for case in json.loads((root/'truth.json').read_text()):
    with open(root/case['name']/'calls/events.tsv') as f:
        events=list(csv.DictReader(f,delimiter='\t'))
    supported=[e for e in events if e.get('status')=='read_supported']
    matched=[e for e in supported if max(abs(int(e[k])-p) for k,p in zip(['a0','a1','b0','b1'],case['coordinates']))<=50]
    passed=bool(matched) if case['expect_supported'] else not supported
    results.append(dict(**case,supported=len(supported),coordinate_matches=len(matched),passed=passed))
(root/'evaluation.json').write_text(json.dumps(results,indent=2)+'\n')
print(f"{sum(r['passed'] for r in results)}/{len(results)} sequence-level smoke cases passed")
if not all(r['passed'] for r in results): sys.exit(1)
