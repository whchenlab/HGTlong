"""Explicit post-call evidence screens; no comparator labels enter this stage."""


def screen_event(event,min_identity=.99,min_coverage=.95,min_reads=3):
    reasons=[]
    if event['status']!='read_supported' or min(int(event['n_'+k]) for k in ('a0','a1','b0','b1'))<min_reads:
        reasons.append('read_support')
    if event.get('qc_status')!='annotation_screen_pass':reasons.append('repeat_or_annotation_qc')
    if event.get('sequence_qc')!='both_blocks_covered':reasons.append('sequence_alignment_qc')
    for key,threshold,reason in [('block_identity',min_identity,'block_identity'),
                                 ('aligned_fraction_A',min_coverage,'coverage_A'),
                                 ('aligned_fraction_B',min_coverage,'coverage_B')]:
        try:value=float(event.get(key,''))
        except (TypeError,ValueError):value=-1
        if not value>=threshold:reasons.append(reason)
    return dict(event,screened_candidate=int(not reasons),screen_fail_reasons=';'.join(reasons))
