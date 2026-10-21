#!/usr/bin/env python3
"""H4: catalogue-wide trait-profile permutation pilot, site-matched treatments."""
import argparse
import csv
import random
import statistics
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--permutations', type=int, default=199)
parser.add_argument('--seed', type=int, default=20261012)
args = parser.parse_args()
assert args.permutations >= 19
rng = random.Random(args.seed)

BASE = Path('results/ecological_strategies')
AB = Path('mapping/competitive_MAGs/MAG_abundance')
OUT = BASE / 'functional_redundancy'
FILES = {
    'Osmotic': 'osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv',
    'Redox': 'redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv',
    'Nitrogen': 'nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv',
    'Carbon_P': 'carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv',
    'Iron': 'iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv',
    'Colonisation': 'colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv',
}
CONTRASTS = {
    'BS': [('NAT','BESE'), ('NAT','BARE'), ('NAT','DEGR'),
           ('BESE','BARE'), ('BESE','DEGR'), ('BARE','DEGR')],
    'RZ': [('NAT','BESE'), ('NAT','BARE'), ('BESE','BARE')],
}

def read_tsv(path):
    with path.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return reader.fieldnames, list(reader)

hd, drows = read_tsv(AB / 'final_373_MAGs_mean_depth_matrix.tsv')
hb, brows = read_tsv(AB / 'final_373_MAGs_breadth_matrix.tsv')
assert hd == hb and hd[0] == 'MAG_ID'
mags = [r['MAG_ID'] for r in drows]
samples = hd[1:]
assert len(mags) == len(set(mags)) == 373
assert len(samples) == 20
assert [r['MAG_ID'] for r in brows] == mags

_, refrows = read_tsv(OUT / 'functional_redundancy_20samples_37capabilities.tsv')
assert len(refrows) == 1480
caps = list(dict.fromkeys((r['GRESP'], r['Capability']) for r in refrows))
assert len(caps) == 37

source = {}
for group, filename in FILES.items():
    _, rows = read_tsv(BASE / filename)
    source[group] = {r['Final_MAG']: r for r in rows}
    assert len(source[group]) == 373 and set(source[group]) == set(mags)

profiles = [
    [j for j, (group, cap) in enumerate(caps)
     if float(source[group][mag][cap]) == 1.0]
    for mag in mags
]

aliases = {'DEG':'DEGR', 'DEGR':'DEGR', 'NAT':'NAT', 'BESE':'BESE', 'BARE':'BARE'}
meta = {}
weights = {}
for sample in samples:
    parts = sample.split('_')
    key = (parts[1], parts[3], aliases[parts[2]])  # site, habitat, treatment
    assert key not in meta.values()
    meta[sample] = key
    for filt, threshold in [('Breadth10', 0.10), ('Breadth50', 0.50)]:
        w = [(i, float(drows[i][sample])) for i in range(373)
             if float(drows[i][sample]) > 0 and float(brows[i][sample]) >= threshold]
        total = sum(v for _, v in w)
        assert total > 0
        weights[(sample, filt)] = [(i, v/total) for i, v in w]

idx = {key: sample for sample, key in meta.items()}
ref = {(r['Sample'], r['Filter'], r['GRESP'], r['Capability']):r for r in refrows}
assert len(ref) == len(refrows)


def compute(assignment):
    """Return (prevalence %, effective carriers, carrier count) for all rows."""
    result = {}
    for key, detected in weights.items():
        wsum = [0.0] * len(caps)
        w2sum = [0.0] * len(caps)
        count = [0] * len(caps)
        for magidx, v in detected:
            for j in profiles[assignment[magidx]]:
                wsum[j] += v
                w2sum[j] += v*v
                count[j] += 1
        result[key] = [(100*a, (a*a/b if b else 0.0), n)
                       for a, b, n in zip(wsum, w2sum, count)]
    return result

observed = compute(list(range(373)))
max_prev_err, max_eff_err = 0., 0.
for (sample, filt), arr in observed.items():
    for j, (group, cap) in enumerate(caps):
        row = ref[(sample, filt, group, cap)]
        p, e, k = arr[j]
        max_prev_err = max(max_prev_err, abs(p-float(row['Weighted_prevalence_pct'])))
        max_eff_err = max(max_eff_err, abs(e-float(row['Effective_carriers'])))
        assert k == int(row['Carrier_richness'])
assert max_prev_err < 0.0001 and max_eff_err < 0.0001, (max_prev_err, max_eff_err)

comparisons = []
for filt in ('Breadth10', 'Breadth50'):
    for habitat, pairs in CONTRASTS.items():
        for first, second in pairs:
            site_pairs = [
                (site, (idx[(site, habitat, first)], filt), (idx[(site, habitat, second)], filt))
                for site in ('ACP', 'BASS', 'TRUG')
                if (site, habitat, first) in idx and (site, habitat, second) in idx
            ]
            assert len(site_pairs) >= 2
            for j, (group, cap) in enumerate(caps):
                comparisons.append((filt, habitat, first, second, group, cap, j, site_pairs))
assert len(comparisons) == 666

def contrasts(values, j, site_pairs):
    site_vals = [(site, tuple(values[k2][j][metric]-values[k1][j][metric] for metric in (0, 1)))
                 for site, k1, k2 in site_pairs]
    means = [statistics.mean(v[metric] for _, v in site_vals) for metric in (0, 1)]
    return means, site_vals

observed_stats = [contrasts(observed, c[6], c[7]) for c in comparisons]
null_stats = [[[], []] for _ in comparisons]
for rep in range(args.permutations):
    assignment = list(range(373))
    rng.shuffle(assignment)  # ONE profile permutation, reused for all 20 samples and two filters
    perm_values = compute(assignment)
    for i, c in enumerate(comparisons):
        means, _ = contrasts(perm_values, c[6], c[7])
        for metric in (0, 1):
            null_stats[i][metric].append(means[metric])
    if (rep+1) % 50 == 0:
        print('Completed randomisations:', rep+1, flush=True)


def quantile(xs, q):
    s = sorted(xs)
    pos = (len(s)-1)*q
    a = int(pos)
    b = min(a+1, len(s)-1)
    return s[a] + (pos-a)*(s[b]-s[a])

records = []
for i, c in enumerate(comparisons):
    filt, habitat, first, second, group, cap, j, site_pairs = c
    means, site_vals = observed_stats[i]
    row = {'Filter':filt, 'Habitat':habitat, 'GRESP':group, 'Capability':cap,
           'Contrast':f'{second}_minus_{first}', 'N_sites':len(site_pairs)}
    for m, tag in enumerate(('Prevalence_pp', 'Effective_carriers')):
        v = means[m]
        null = null_stats[i][m]
        mu = statistics.mean(null)
        # empirical two-sided extremeness around the null distribution's centre
        p = (1+sum(abs(x-mu) >= abs(v-mu)-1e-12 for x in null))/(len(null)+1)
        row.update({
            f'Observed_delta_{tag}':round(v, 6),
            f'Null_mean_delta_{tag}':round(mu, 6),
            f'Null_95_low_{tag}':round(quantile(null, 0.025), 6),
            f'Null_95_high_{tag}':round(quantile(null, 0.975), 6),
            f'Empirical_p_{tag}':round(p, 6),
            f'Positive_sites_{tag}':sum(x[m] > 1e-12 for _,x in site_vals),
            f'Negative_sites_{tag}':sum(x[m] < -1e-12 for _,x in site_vals),
        })
        for site in ('ACP', 'BASS', 'TRUG'):
            vv = next((x[m] for s, x in site_vals if s == site), None)
            row[f'{site}_delta_{tag}'] = '' if vv is None else round(vv, 6)
    row['Sites_both_2plus_carriers'] = sum(
        observed[k1][j][2] >= 2 and observed[k2][j][2] >= 2
        for _, k1, k2 in site_pairs
    )
    records.append(row)

OUT.mkdir(parents=True, exist_ok=True)
out = OUT / f'functional_catalogue_null_{args.permutations}perms.tsv'
with out.open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(records)
print('SUCCESS: Catalogue-wide trait allocation null')
print('MAGs:', len(mags), '| samples:', len(samples), '| capabilities:', len(caps))
print('Randomisations:', args.permutations, '| matched contrast rows:', len(records))
print('Maximum observed-data QC error (prevalence / effective):', max_prev_err, max_eff_err)
print('Output:', out)
print('NOTE: Empirical p-values test trait allocation, not biological treatment significance.')

