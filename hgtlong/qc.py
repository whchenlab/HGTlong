"""Post-call annotation overlaps. This stage cannot create candidate events."""
import csv
from collections import defaultdict
from .core import union_length


def annotate_masks(events, bed_path, screened_contigs):
    masks=defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            fields=line.rstrip().split('\t')
            if len(fields)<3:raise ValueError('Mask BED needs at least three columns')
            contig,start,end=fields[0],int(fields[1]),int(fields[2])
            if start<0 or end<=start:raise ValueError('Invalid mask interval')
            masks[contig].append((start,end,fields[3] if len(fields)>3 else 'repeat'))
    mask_classes=','.join(sorted({label for values in masks.values() for _,_,label in values})) or 'unspecified'
    for event in events:
        labels=set()
        for side in ['A','B']:
            lo,hi=int(event[side.lower()+'0']),int(event[side.lower()+'1'])
            overlaps=[]
            for start,end,label in masks[event['ctg'+side]]:
                if min(end,hi)>max(start,lo):
                    overlaps.append((max(start,lo),min(end,hi)));labels.add(label)
            event['mask_bp_'+side]=union_length(overlaps)
        covered=all(event['ctg'+s] in screened_contigs for s in ['A','B'])
        event['repeat_screen']=','.join(sorted(labels)) if labels else 'no_overlap_in_supplied_masks' if covered else 'incomplete_annotation'
        event['qc_status']='repeat_associated' if labels else 'annotation_screen_pass' if covered else 'not_fully_assessed'
        event['mask_classes']=mask_classes
        yield event
