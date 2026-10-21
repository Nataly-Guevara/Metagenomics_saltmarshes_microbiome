#!/usr/bin/env python3
"""Compare whole-MAG and intact-profile capability-count-weighted shared abundance.

Input matrices and carrier calls must be prepared by the existing Gadi workflow.
One global intact-profile permutation per iteration, shared across samples/filters.
"""
import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

AB = Path('mapping/competitive_MAGs/MAG_abundance')
BASE = Path('results/ecological_strategies/functional_carrier_taxonomy')
OUT = BASE / 'carrier_turnover'
QC = Path('results/ecological_strategies/final_validation/HQ_37capability_catalogue.tsv')
NPERM = 999
SEED = 20261014
CONTRASTS = [('NAT','BESE'), ('NAT','BARE'), ('BESE','BARE')]
SITES = ['ACP', 'BASS', 'TRUG']
HABITATS = ['BS','RZ']
FILTERS = [('Breadth10',0.10), ('Breadth50',0.50)]

def read_tsv(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(p, rows):
    assert rows
    with p.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(rows)

def percentile(a,p):
    s=sorted(a)
    x=(len(s)-1)*p
    lo=int(math.floor(x)); hi=int(math.ceil(x))
    return s[lo]*(hi-x)+s[hi]*(x-lo) if hi>lo else s[lo]

def bh(p_values):
    n=len(p_values)
    order=sorted(range(n),key=lambda x:p_values[x])
    result=[0.0]*n; cumulative=1.0
    for rank in range(n-1,-1,-1):
        i=order[rank]
        cumulative=min(cumulative,p_values[i]*n/(rank+1))
        result[i]=cumulative
    return result

D=read_tsv(AB/'final_373_MAGs_mean_depth_matrix.tsv')
B=read_tsv(AB/'final_373_MAGs_breadth_matrix.tsv')
carrier=read_tsv(BASE/'MAG_capability_carrier_taxonomy.tsv')
catalogue=read_tsv(QC)
assert len(D)==373 and len(B)==373
ids=[r['MAG_ID'] for r in D]
assert len(set(ids))==373
by_breadth={r['MAG_ID']:r for r in B}
assert set(ids)==set(by_breadth)
samples=[x for x in D[0] if x!='MAG_ID']
assert len(samples)==20 and set(B[0])==set(D[0])
capabilities={(r['GRESP'],r['Capability']) for r in catalogue}
assert len(capabilities)==37
profiles={m:set() for m in ids}
for r in carrier:
    mag=r['MAG_ID']; c=(r['GRESP'],r['Capability'])
    assert mag in profiles and c in capabilities
    if c in profiles[mag]:
        raise ValueError('Duplicate genome-capability record: '+mag+' '+str(c))
    profiles[mag].add(c)
# Each permutation shuffles entire intact genome profiles; only profile sizes
# affect THIS statistic, but correlated capability calls remain together.
original_k={mag:len(profiles[mag]) for mag in ids}
profile_list=[frozenset(profiles[m]) for m in ids]

cells={}
for s in samples:
    fields=s.split('_')
    assert len(fields)>=4
    key=tuple(fields[1:4])
    assert key not in cells
    cells[key]=s

all_weights={}
for filter_name,threshold in FILTERS:
    weights={}
    for s in samples:
        retained={
            r['MAG_ID']: float(r[s]) for r in D
            if float(r[s])>0 and float(by_breadth[r['MAG_ID']][s])>=threshold
        }
        total=sum(retained.values())
        assert total>0 and all(x>=0 for x in retained.values())
        weights[s]={mag:v/total for mag,v in retained.items()}
    all_weights[filter_name]=weights

def one_direction(w, shared, k):
    W=sum(w[m] for m in shared)  # 0..1, includes zero-capability MAGs
    weighted_total=sum(v*k[m] for m,v in w.items())
    if weighted_total<=0:
        raise ValueError('All detected MAGs have zero curated capabilities')
    weighted_shared=sum(w[m]*k[m] for m in shared)
    C=weighted_shared/weighted_total
    mean_shared=weighted_shared/W if W>0 else float('nan')
    mean_nonshared=(weighted_total-weighted_shared)/(1-W) if W<1 else float('nan')
    return W*100,C*100,(C-W)*100,mean_shared,mean_nonshared

pairs=[]
for name,_threshold in FILTERS:
    weights=all_weights[name]
    for A,Bname in CONTRASTS:
        for site in SITES:
            for habitat in HABITATS:
                ca=(site,A,habitat); cb=(site,Bname,habitat)
                if ca not in cells or cb not in cells:
                    continue
                wa=weights[cells[ca]]; wb=weights[cells[cb]]
                shared=set(wa)&set(wb)
                A_vals=one_direction(wa,shared,original_k)
                B_vals=one_direction(wb,shared,original_k)
                pairs.append({
                    'Filter':name,'Site':site,'Habitat':habitat,
                    'Contrast':A+'_vs_'+Bname,
                    'A':A,'B':Bname,
                    'wa':wa,'wb':wb,'shared':shared,
                    'A_values':A_vals,'B_values':B_vals,
                    'Observed_pair_C_minus_W_pp':mean([A_vals[2],B_vals[2]]),
                    'null_A':[],'null_B':[],'null_pair':[]
                })
assert len(pairs)==32, len(pairs)

# Cross-check independently calculated whole-MAG shared abundance.
reference=OUT/'Whole_MAG_vs_capability_shared_abundance.tsv'
if reference.exists():
    baseline={(r['Filter'],r['Site'],r['Habitat'],r['Contrast']):r
              for r in read_tsv(reference)}
    for r in pairs:
        key=tuple(r[x] for x in ('Filter','Site','Habitat','Contrast'))
        check=baseline[key]
        assert abs(r['A_values'][0]-float(check['Whole_MAG_shared_A_pct']))<1e-7
        assert abs(r['B_values'][0]-float(check['Whole_MAG_shared_B_pct']))<1e-7
    print('Whole-MAG shared abundance verified against existing baseline for all 32 pairs.')

rng=random.Random(SEED)
for iteration in range(NPERM):
    shuffled=profile_list.copy()
    rng.shuffle(shuffled)  # same genome-profile assignment in all samples
    perm_k={m:len(p) for m,p in zip(ids,shuffled)}
    for r in pairs:
        a=one_direction(r['wa'],r['shared'],perm_k)[2]
        b=one_direction(r['wb'],r['shared'],perm_k)[2]
        r['null_A'].append(a)
        r['null_B'].append(b)
        r['null_pair'].append((a+b)/2)

rows=[]
for r in pairs:
    a=r['A_values']; b=r['B_values']
    null=r['null_pair']; expected=mean(null)
    obs=r['Observed_pair_C_minus_W_pp']
    two=(1+sum(abs(x-expected)>=abs(obs-expected)-1e-12 for x in null))/(NPERM+1)
    upper=(1+sum(x>=obs-1e-12 for x in null))/(NPERM+1)
    rows.append({
        'Filter':r['Filter'],'Site':r['Site'],'Habitat':r['Habitat'],
        'Contrast':r['Contrast'],'Treatment_A':r['A'],'Treatment_B':r['B'],
        'N_shared_MAGs':len(r['shared']),
        'Whole_shared_A_W_pct':a[0], 'Capability_shared_A_C_pct':a[1],
        'C_minus_W_A_pp':a[2], 'Shared_mean_k_A':a[3],
        'Nonshared_mean_k_A':a[4], 'Null_expected_A_pp':mean(r['null_A']),
        'Whole_shared_B_W_pct':b[0], 'Capability_shared_B_C_pct':b[1],
        'C_minus_W_B_pp':b[2], 'Shared_mean_k_B':b[3],
        'Nonshared_mean_k_B':b[4], 'Null_expected_B_pp':mean(r['null_B']),
        'Pair_W_pct':mean([a[0],b[0]]),
        'Pair_C_pct':mean([a[1],b[1]]),
        'Pair_C_minus_W_pp':obs,
        'Null_pair_expected_pp':expected,
        'Null_pair_2.5pct_pp':percentile(null,0.025),
        'Null_pair_97.5pct_pp':percentile(null,0.975),
        'Observed_minus_null_pp':obs-expected,
        'One_sided_upper_p':upper,
        'Two_sided_centered_p':two,
        'N_permutations':NPERM
    })

# BH correction within each breadth filter across all site/contrast/habitat
# comparisons; breadth50 is a sensitivity analysis, not independent evidence.
for filt,_ in FILTERS:
    positions=[i for i,r in enumerate(rows) if r['Filter']==filt]
    q=bh([rows[i]['Two_sided_centered_p'] for i in positions])
    for i,qi in zip(positions,q):rows[i]['BH_q_within_breadth']=qi
write_tsv(OUT/'Shared_MAG_capability_count_null_pair_results.tsv',rows)

summ=[]
for filt,_ in FILTERS:
    for A,Bname in CONTRASTS:
        contrast=A+'_vs_'+Bname
        for habitat in HABITATS:
            sub=[r for r in rows if r['Filter']==filt and r['Contrast']==contrast and r['Habitat']==habitat]
            if not sub: continue
            summ.append({
                'Filter':filt,'Contrast':contrast,'Habitat':habitat,
                'N_sites':len(sub),'Sites':','.join(r['Site'] for r in sub),
                'Mean_W_pct':mean(r['Pair_W_pct'] for r in sub),
                'Mean_C_pct':mean(r['Pair_C_pct'] for r in sub),
                'Mean_C_minus_W_pp':mean(r['Pair_C_minus_W_pp'] for r in sub),
                'Mean_null_expected_pp':mean(r['Null_pair_expected_pp'] for r in sub),
                'Mean_observed_minus_null_pp':mean(r['Observed_minus_null_pp'] for r in sub)
            })
write_tsv(OUT/'Shared_MAG_capability_count_null_summary.tsv',summ)

print('\n=== PRIMARY 10% BREADTH: SITE-SPECIFIC C-W ===')
print('Site Habitat Contrast    W%     C%    C-Wpp   Null Epp    centered p   BH q')
for r in rows:
    if r['Filter']!='Breadth10': continue
    print(f"{r['Site']:<5}{r['Habitat']:<8}{r['Contrast']:<14}"
          f"{r['Pair_W_pct']:6.1f} {r['Pair_C_pct']:6.1f} "
          f"{r['Pair_C_minus_W_pp']:+7.2f} {r['Null_pair_expected_pp']:+8.2f} "
          f"{r['Two_sided_centered_p']:7.3f} {r['BH_q_within_breadth']:7.3f}")
print('\n=== SUMMARY BY COMPARISON AND HABITAT ===')
for s in summ:
    print(f"{s['Filter']} {s['Contrast']} {s['Habitat']} "
          f"sites={s['N_sites']}: W={s['Mean_W_pct']:.1f}% "
          f"C={s['Mean_C_pct']:.1f}% "
          f"C-W={s['Mean_C_minus_W_pp']:+.2f}pp "
          f"null E={s['Mean_null_expected_pp']:+.2f}pp")
print('\nSUCCESS: Intact-profile null,',NPERM,'permutations; seed',SEED)
print('Output:',OUT/'Shared_MAG_capability_count_null_pair_results.tsv')
print('Output:',OUT/'Shared_MAG_capability_count_null_summary.tsv')
print('Interpretation: null allocation ignores phylogenetic and genome-quality structure.')
print('Non-independent sample/contrast results are descriptive, not habitat-effect tests.')
