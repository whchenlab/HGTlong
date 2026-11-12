"""Pairwise comparator geometry; HDMI is not a truth label."""

def canonical_geometry(row):
    a=(row['ctgA'],int(row['a0']),int(row['a1']))
    b=(row['ctgB'],int(row['b0']),int(row['b1']))
    return tuple(sorted([a,b]))+(row['orientation'],)


def compare_geometry(call,reference):
    if {call['ctgA'],call['ctgB']}!={reference['ctgA'],reference['ctgB']}:
        return None
    keys=['a0','a1','b0','b1'] if call['ctgA']==reference['ctgA'] else ['b0','b1','a0','a1']
    bounds=[int(reference[k]) for k in keys]
    observed=[int(call[k]) for k in ['a0','a1','b0','b1']]
    overlap=[max(0,min(observed[i+1],bounds[i+1])-max(observed[i],bounds[i])) for i in [0,2]]
    orientation=call['orientation']==reference['orientation']
    distance=max(abs(a-b) for a,b in zip(observed,bounds))
    return dict(exact=orientation and distance==0,bilateral=orientation and all(overlap),orientation_match=orientation,max_boundary_distance=distance,overlap_A_bp=overlap[0],overlap_B_bp=overlap[1],fraction_call_A=overlap[0]/(observed[1]-observed[0]),fraction_call_B=overlap[1]/(observed[3]-observed[2]),fraction_HDMI_A=overlap[0]/(bounds[1]-bounds[0]),fraction_HDMI_B=overlap[1]/(bounds[3]-bounds[2]))


def overlap_category(matches):
    if any(m['final_retained'] and m['exact'] for m in matches):return 'identical_to_HDMI_final'
    if any(m['final_retained'] and m['bilateral'] for m in matches):return 'partial_overlap_with_HDMI_final'
    if any(m['bilateral'] for m in matches):return 'overlap_HDMI_raw_not_final'
    return 'no_bilateral_HDMI_raw_overlap'
