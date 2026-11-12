"""Coordinate helpers for exploratory annotation-context analyses."""

def merge_intervals(intervals):
    merged=[]
    for start,end in sorted(intervals):
        if start>=end:raise ValueError('Empty or reversed interval')
        if merged and start<=merged[-1][1]:merged[-1]=(merged[-1][0],max(end,merged[-1][1]))
        else:merged.append((start,end))
    return merged


def shift_intervals(intervals,shift,length):
    pieces=[]
    for start,end in intervals:
        size=end-start
        if not 0<=start<end<=length:raise ValueError('Interval outside contig')
        start=(start+shift)%length;end=start+size
        if end<=length:pieces.append((start,end))
        else:pieces.extend([(start,length),(0,end-length)])
    return merge_intervals(pieces)


def bh_adjust(values):
    order=sorted(range(len(values)),key=lambda i:values[i]);result=[1.]*len(values);running=1.
    for j in range(len(order)-1,-1,-1):
        i=order[j];running=min(running,values[i]*len(values)/(j+1));result[i]=min(1.,running)
    return result


def interval_distance(a,b):
    return max(0,a[0]-b[1],b[0]-a[1])


def boundary_clusters(points,width=1000):
    """Fixed maximum diameter, not chained proximity. Points must have position."""
    groups=[]
    for point in sorted(points,key=lambda p:p['position']):
        if not groups or point['position']-groups[-1][0]['position']>width:groups.append([point])
        else:groups[-1].append(point)
    return groups


def parse_cog_categories(lines):
    """NCBI COG2024 permits category rows without a display-color column."""
    categories={}
    for line in lines:
        fields=line.rstrip('\r\n').split('\t')
        if len(fields)>=3 and len(fields[0])==1 and fields[0].isalpha():
            categories[fields[0]]=dict(name=fields[-1],color='#'+fields[2].lstrip('#') if len(fields)>=4 else '#c78bbd')
    return categories
