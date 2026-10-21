#!/usr/bin/env python3
"""Taxon abundance x carrier fraction decomposition for saltmarsh MAGs.

No third-party packages. Run from the Gadi metagenomics project directory.
All 37 capabilities: phylum profiles/decomposition.
Seven previously examined focal capabilities: order and individual-MAG resolution.
Both domains, both breadth filters, all-MAG and >=90%/<5% subsets.
This is exact descriptive accounting, not a causal or significance test.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime
from pathlib import Path
from statistics import mean

GROUPS = {
    "Osmotic": ("osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
        "Capability_Trehalose Capability_Compatible_solute_uptake Capability_Ectoine_biosynthesis Capability_Mechanosensitive_protection Capability_Ion_homeostasis"),
    "Redox": ("redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
        "Capability_Aerobic_respiration Capability_LowO2_respiration Capability_Nitrate_respiration Capability_Nitrite_to_NO Capability_DNRA Capability_NO_reduction Capability_N2O_reduction Capability_Dissimilatory_sulfate_reduction_core Capability_Reduced_sulfur_oxidation Capability_Sulfite_oxidation"),
    "Nitrogen": ("nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
        "Capability_Nitrogen_fixation Capability_Ammonium_uptake Capability_Assimilatory_nitrate_reduction Capability_Urea_utilisation Capability_DNRA Capability_Complete_denitrification Capability_Nitrification Capability_Anammox"),
    "Carbon_P": ("carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
        "Capability_Carbon_fixation Capability_C1_utilisation Capability_High_affinity_phosphate_uptake Capability_Organic_P_mobilisation Capability_Phosphonate_acquisition Capability_PolyP_turnover"),
    "Iron": ("iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
        "Ferrous_iron_uptake Direct_ferric_iron_uptake Siderophore_biosynthesis Iron_complex_uptake"),
    "Colonisation": ("colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv",
        "Chemotaxis Motility Attachment_biofilm_matrix Direct_plant_interaction")
}
FOCAL = [
    ("Osmotic", "Capability_Trehalose"),
    ("Carbon_P", "Capability_PolyP_turnover"),
    ("Nitrogen", "Capability_Ammonium_uptake"),
    ("Iron", "Iron_complex_uptake"),
    ("Redox", "Capability_LowO2_respiration"),
    ("Redox", "Capability_Nitrate_respiration"),
    ("Redox", "Capability_N2O_reduction")
]
# A is subtracted from B. Explicit names prevent sign ambiguities.
CONTRASTS = [("NAT", "BESE"), ("NAT", "BARE"), ("NAT", "DEGR"),
             ("BARE", "BESE"), ("DEGR", "BESE"), ("DEGR", "BARE")]
SITES = ("ACP", "BASS", "TRUG")
HABITATS = ("BS", "RZ")
CONDITIONS = ("NAT", "BESE", "BARE", "DEGR")
FILTERS = (("Breadth10", 0.10), ("Breadth50", 0.50))
BASE = Path("results/ecological_strategies")
AB = Path("mapping/competitive_MAGs/MAG_abundance")
IDFIELDS = ["Filter", "Subset", "Comparison_type", "Stratum", "Site", "Contrast",
            "Sample_A", "Sample_B", "Treatment_A", "Treatment_B", "Habitat_A", "Habitat_B"]
TAXFIELDS = ["Rank", "Taxon_ID", "Domain", "Phylum", "Order"]
COMPONENTS = ["Representation_pp", "Within_lineage_mixture_pp", "One_sided_lineage_pp"]


def read_tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        raise ValueError("Empty file: " + str(path))
    return rows, fields


def keyed(rows, field, path):
    if field not in rows[0]:
        raise ValueError("Missing %s in %s" % (field, path))
    result = {}
    for r in rows:
        k = r[field]
        if not k or k in result:
            raise ValueError("Missing/duplicate identifier in %s: %s" % (path, k))
        result[k] = r
    return result


def numeric(value, label, lo=0.0, hi=None):
    x = float(value)
    if not math.isfinite(x) or x < lo or (hi is not None and x > hi):
        raise ValueError("Invalid numeric value %r: %s" % (value, label))
    return x


def binary(value):
    s = str(value).strip().lower()
    if s in ("1", "1.0", "true"): return 1
    if s in ("0", "0.0", "false"): return 0
    raise ValueError("Invalid binary capability call: " + str(value))


def close(a, b, label, tolerance=2e-8):
    if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > tolerance:
        raise ValueError("Validation failed: %s; %r versus %r" % (label, a, b))


def writer(stack, path, fields):
    opener = gzip.open if path.name.endswith(".gz") else open
    f = stack.enter_context(opener(path, "wt", newline="", encoding="utf-8"))
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="raise")
    w.writeheader()
    return w


def taxonomy_record(r):
    terms = {}
    for token in r["GTDB_classification"].split(";"):
        token = token.strip()
        if "__" in token:
            k, v = token.split("__", 1)
            terms[k] = v
    domain = terms.get("d", "") or r.get("GTDB_domain", "")
    domain = domain.replace("d__", "")
    if domain not in ("Bacteria", "Archaea"):
        raise ValueError("Unexpected domain: " + domain)
    if r.get("GTDB_domain") and r["GTDB_domain"].replace("d__", "") != domain:
        raise ValueError("Domain disagrees with full taxonomy")
    terms["d"] = domain
    # Full lineage paths prevent collisions among unnamed/order labels.
    pkey = ";".join(k + "__" + terms.get(k, "") for k in ("d", "p"))
    okey = ";".join(k + "__" + terms.get(k, "") for k in ("d", "p", "c", "o"))
    return {"Domain": domain, "Phylum": terms.get("p", "") or "Unclassified phylum",
            "Order": terms.get("o", "") or "Unclassified order", "Phylum_ID": pkey,
            "Order_ID": okey, "GTDB_classification": r["GTDB_classification"]}


def decompose(pa, fa, pb, fb):
    """P in percent, F in [0,1] or None for an undetected lineage.

    Return midpoint decomposition where both P>0. One-sided lineages
    are reported separately rather than inventing an F for the absent side.
    """
    ya = pa * fa if pa > 0 else 0.0
    yb = pb * fb if pb > 0 else 0.0
    if pa > 0 and pb > 0:
        representation = (pb - pa) * (fa + fb) / 2.0
        mixture = (fb - fa) * (pa + pb) / 2.0
        single = 0.0
        status = "Detected_both"
    else:
        representation = mixture = 0.0
        single = yb - ya
        status = "Detected_B_only" if pb > 0 else "Detected_A_only" if pa > 0 else "Detected_neither"
    error = representation + mixture + single - (yb - ya)
    close(error, 0.0, "product decomposition")
    return representation, mixture, single, status, error


def self_test():
    # Increased lineage abundance at fixed carrier fraction: +5 pp.
    x = decompose(20, .5, 30, .5); close(x[0], 5, "representation test"); close(x[1], 0, "mixture test")
    # Increased carrier fraction at fixed lineage abundance: +4 pp.
    x = decompose(20, .5, 20, .7); close(x[0], 0, "representation test"); close(x[1], 4, "mixture test")
    # Exact symmetric allocation when both quantities change.
    x = decompose(20, .5, 30, .7); close(x[0], 6, "midpoint R"); close(x[1], 5, "midpoint M")
    y = decompose(30, .7, 20, .5)
    for i in range(3): close(x[i], -y[i], "direction reversal")
    x = decompose(0, None, 30, .7); close(x[2], 21, "one-sided emergence")
    x = decompose(20, 0.0, 20, 0.0); close(sum(x[:3]), 0, "noncarrier lineage")
    x = decompose(0, None, 0, None); close(sum(x[:3]), 0, "both absent")
    print("SELF-TEST PASSED: products, midpoint identity, sign reversal, zero carriers and absent lineages.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path.cwd(), help="Gadi project root")
    ap.add_argument("--out", type=Path, default=None, help="New output directory; refuses to overwrite")
    ap.add_argument("--self-test", action="store_true", help="Test accounting functions without input data")
    args = ap.parse_args()
    if args.self_test:
        self_test(); return
    root = args.root.expanduser().resolve()
    sources = []
    def load(relative, idcol=None):
        p = root / relative
        rows, fields = read_tsv(p)
        sources.append(p)
        return (keyed(rows, idcol, p) if idcol else rows), fields
    depth_rows, dh = load(AB / "final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID")
    breadth_rows, bh = load(AB / "final_373_MAGs_breadth_matrix.tsv", "MAG_ID")
    taxrows, _ = load(Path("results/final_summaries/combined_binners/final_373_MAGs_GTDB_taxonomy.tsv"), "Final_MAG")
    mags = sorted(depth_rows)
    samples = [x for x in dh if x != "MAG_ID"]
    if len(mags) != 373 or len(samples) != 20:
        raise ValueError("Expected 373 MAGs and 20 metagenomes")
    if set(mags) != set(breadth_rows) or set(mags) != set(taxrows) or set(dh) != set(bh):
        raise ValueError("Depth, breadth and taxonomy identifiers/columns do not match")
    taxonomy = {m: taxonomy_record(taxrows[m]) for m in mags}
    domains = Counter(t["Domain"] for t in taxonomy.values())
    if domains != Counter(Bacteria=352, Archaea=21):
        raise ValueError("Unexpected domain counts: " + str(domains))
    depth = {s: {m: numeric(depth_rows[m][s], "depth") for m in mags} for s in samples}
    breadth = {s: {m: numeric(breadth_rows[m][s], "breadth", hi=1) for m in mags} for s in samples}
    cells, metadata = {}, {}
    for s in samples:
        p = s.split("_")
        if len(p) < 4: raise ValueError("Unrecognised sample ID: " + s)
        site, condition, habitat = p[1:4]
        if site not in SITES or condition not in CONDITIONS or habitat not in HABITATS:
            raise ValueError("Unrecognised metadata: " + s)
        key = (site, condition, habitat)
        if key in cells: raise ValueError("Duplicate site x condition x habitat")
        cells[key] = s
        metadata[s] = {"Site": site, "Treatment": condition, "Habitat": habitat}
    expected = {(s,t,h) for s in SITES for t in CONDITIONS for h in HABITATS
                if not (h == "RZ" and t == "DEGR") and (s,t,h) != ("TRUG","BARE","BS")}
    if set(cells) != expected:
        raise ValueError("Sample design differs from documented 20 metagenomes")
    calls, quality = {}, None
    for group, (name, fields) in GROUPS.items():
        rows, _ = load(BASE / name, "Final_MAG")
        if set(rows) != set(mags): raise ValueError("MAG IDs differ in " + group)
        for field in fields.split():
            calls[(group, field)] = {m: binary(rows[m][field]) for m in mags}
        if group == "Iron":
            quality = {m: (numeric(rows[m]["Completeness"], "completeness", hi=100),
                           numeric(rows[m]["Contamination"], "contamination")) for m in mags}
    if len(calls) != 37 or not set(FOCAL).issubset(calls):
        raise ValueError("Capability set inconsistent")
    features = list(calls)
    carriers = {c: [m for m in mags if calls[c][m]] for c in features}
    hq = {m for m in mags if quality[m][0] >= 90 and quality[m][1] < 5}
    if not hq: raise ValueError("No genomes meet requested quality thresholds")
    print("Catalogue:", len(mags), dict(domains))
    print("Strict quality subset:", len(hq), dict(Counter(taxonomy[m]["Domain"] for m in hq)), flush=True)
    # Optional direct consistency check against previously exported carrier IDs.
    exported_carriers = root / BASE / "functional_carrier_taxonomy/MAG_capability_carrier_taxonomy.tsv"
    checks = []
    if exported_carriers.exists():
        erows, _ = read_tsv(exported_carriers); sources.append(exported_carriers)
        actual = {(r["GRESP"], r["Capability"], r["MAG_ID"]) for r in erows}
        target = {(g,c,m) for (g,c), ms in carriers.items() for m in ms}
        if actual != target or len(actual) != len(erows):
            raise ValueError("Final capability calls differ from carrier export; reconcile versions before analysis")
        checks.append("Final binary calls agree with existing carrier-ID export.")
    # Named groups include domain and full GTDB path, keeping archaeal carriers.
    members, taxinfo = {}, {}
    for rank in ("Phylum", "Order"):
        groups = defaultdict(list)
        for m in mags: groups[taxonomy[m][rank + "_ID"]].append(m)
        members[rank] = dict(groups)
        for key, ms in groups.items():
            t = taxonomy[ms[0]]
            taxinfo[(rank,key)] = {"Rank": rank, "Taxon_ID": key, "Domain": t["Domain"],
                                   "Phylum": t["Phylum"], "Order": t["Order"] if rank == "Order" else ""}
    # All retained MAGs, including noncarriers, remain in abundance denominators.
    weights, coverage, profile, prevalence = {}, [], {}, {}
    for filt, cutoff in FILTERS:
        for s in samples:
            kept_all = {m: depth[s][m] for m in mags if depth[s][m] > 0 and breadth[s][m] >= cutoff}
            total_all = sum(kept_all.values())
            if total_all <= 0: raise ValueError("No detected full-catalogue MAGs: " + s)
            for subset, eligible in (("All373", set(mags)), ("HQ90_C5", hq)):
                kept = {m:d for m,d in kept_all.items() if m in eligible}
                total = sum(kept.values())
                w = {m:d/total for m,d in kept.items()} if total else {}
                weights[(filt,subset,s)] = w
                rec = dict(Filter=filt, Subset=subset, Sample=s, **metadata[s], Catalogue_MAGs=len(eligible),
                           Detected_MAGs=len(w), Detected_bacterial_MAGs=sum(taxonomy[m]["Domain"]=="Bacteria" for m in w),
                           Detected_archaeal_MAGs=sum(taxonomy[m]["Domain"]=="Archaea" for m in w),
                           Retained_total_depth=total, All373_detected_total_depth=total_all,
                           Retained_share_of_All373_depth_pct=100*total/total_all,
                           Archaeal_abundance_pct=100*sum(v for m,v in w.items() if taxonomy[m]["Domain"]=="Archaea") if total else "",
                           Status="Evaluable" if total else "No_retained_MAGs")
                coverage.append(rec)
                for c in features:
                    prevalence[(filt,subset,s,c)] = 100*sum(w.get(m,0) for m in carriers[c]) if total else None
                if not total: continue
                for rank in members:
                    fset = features if rank == "Phylum" else FOCAL
                    for key, ms in members[rank].items():
                        detected = [m for m in ms if m in w]
                        p = 100*sum(w[m] for m in detected)
                        for c in fset:
                            cm = [m for m in detected if calls[c][m]]
                            y = 100*sum(w[m] for m in cm)
                            f = y/p if p > 0 else None
                            profile[(filt,subset,s,rank,key,c)] = (p,f,y,len(detected),len(cm))
                    for c in fset:
                        reconstructed = sum(profile[(filt,subset,s,rank,t,c)][2] for t in members[rank])
                        close(reconstructed, prevalence[(filt,subset,s,c)], "taxa sum to capability prevalence")
    # Compare independently reconstructed prevalence with existing exported values.
    refs = [(BASE / "community_weighted_GRESP/community_capabilities_20_samples_long.tsv", "Weighted_prevalence_pct", False),
            (BASE / "final_validation/HQ_37capability_sample_prevalence.tsv", "Prevalence_pct", True)]
    for rel, col, has_subset in refs:
        p = root / rel
        if not p.exists():
            checks.append("Not available; optional baseline check skipped: " + str(rel)); continue
        rr, _ = read_tsv(p); sources.append(p); count = 0; err = 0
        for r in rr:
            key = (r["Filter"], r["Subset"] if has_subset else "All373", r["Sample"], (r["GRESP"],r["Capability"]))
            if key not in prevalence or prevalence[key] is None: continue
            v = numeric(r[col], "reference prevalence", hi=100.000001)
            diff = abs(v - prevalence[key]); err = max(err,diff)
            close(v, prevalence[key], "existing prevalence " + str(key), tolerance=2e-5); count += 1
        checks.append("%s: %d prevalence values verified; max absolute difference %.3g pp." % (rel.name,count,err))
    output = (args.out if args.out else root / BASE / "lineage_capability_decomposition" /
              ("run_" + datetime.now().strftime("%Y%m%d_%H%M%S")))
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    # Matched comparisons never impute missing combinations.
    comparisons = []
    for site in SITES:
        for hab in HABITATS:
            for a,b in CONTRASTS:
                if (site,a,hab) in cells and (site,b,hab) in cells:
                    comparisons.append(("Condition",hab,site,b+"_minus_"+a,cells[(site,a,hab)],cells[(site,b,hab)]))
        for condition in ("NAT","BESE","BARE"):
            if (site,condition,"BS") in cells and (site,condition,"RZ") in cells:
                comparisons.append(("Habitat",condition,site,"RZ_minus_BS",cells[(site,condition,"BS")],cells[(site,condition,"RZ")]))
    if len(comparisons) != 32: raise ValueError("Unexpected number of matched comparisons")
    site_rows, top_rows, skipped, max_error = [], [], [], 0.0
    focal_taxon_groups = defaultdict(list)
    with ExitStack() as stack:
        cw = writer(stack, output/"Sample_catalogue_coverage.tsv", list(coverage[0])); cw.writerows(coverage)
        pw = writer(stack, output/"Phylum_capability_profiles.tsv.gz", ["Filter","Subset","Sample","Site","Treatment","Habitat"]+TAXFIELDS+
                    ["GRESP","Capability","Taxon_abundance_pct","Within_taxon_carrier_fraction","Within_taxon_carrier_pct",
                     "Capability_contribution_pp","Community_capability_prevalence_pct","Taxon_share_of_capability_pct","Detected_taxon_MAGs","Detected_carrier_MAGs"])
        ow = writer(stack, output/"Focal_order_capability_profiles.tsv.gz", pw.fieldnames)
        for (filt,subset,s,rank,t,c), (p,f,y,nt,nc) in profile.items():
            q = prevalence[(filt,subset,s,c)]
            row = dict(Filter=filt,Subset=subset,Sample=s,**metadata[s],**taxinfo[(rank,t)],GRESP=c[0],Capability=c[1],
                       Taxon_abundance_pct=p,Within_taxon_carrier_fraction=f if f is not None else "",
                       Within_taxon_carrier_pct=100*f if f is not None else "",Capability_contribution_pp=y,
                       Community_capability_prevalence_pct=q,Taxon_share_of_capability_pct=100*y/q if q else "",
                       Detected_taxon_MAGs=nt,Detected_carrier_MAGs=nc)
            (pw if rank=="Phylum" else ow).writerow(row)
        mw = writer(stack, output/"Focal_MAG_sample_contributions.tsv.gz", ["Filter","Subset","Sample","Site","Treatment","Habitat",
                   "MAG_ID","Domain","Phylum","Order","GTDB_classification","Completeness","Contamination","GRESP","Capability",
                   "Mean_depth","Genome_breadth","Detected","Community_abundance_pct","Capability_contribution_pp"])
        for (filt,subset,s), w in weights.items():
            eligible = set(mags) if subset=="All373" else hq
            for c in FOCAL:
                for m in carriers[c]:
                    if m not in eligible: continue
                    t = taxonomy[m]
                    mw.writerow(dict(Filter=filt,Subset=subset,Sample=s,**metadata[s],MAG_ID=m,
                        **{k:t[k] for k in ("Domain","Phylum","Order","GTDB_classification")},
                        Completeness=quality[m][0],Contamination=quality[m][1],GRESP=c[0],Capability=c[1],
                        Mean_depth=depth[s][m],Genome_breadth=breadth[s][m],Detected=int(m in w),
                        Community_abundance_pct=100*w.get(m,0) if w else "",Capability_contribution_pp=100*w.get(m,0) if w else ""))
        tw = writer(stack, output/"Matched_taxon_decomposition.tsv.gz", IDFIELDS+TAXFIELDS+["GRESP","Capability",
                    "Taxon_abundance_A_pct","Taxon_abundance_B_pct","Within_taxon_carrier_A_fraction","Within_taxon_carrier_B_fraction",
                    "Taxon_contribution_A_pp","Taxon_contribution_B_pp","Net_taxon_difference_pp"]+COMPONENTS+
                    ["Taxon_detection_status","Detected_MAGs_A","Detected_MAGs_B","Carrier_MAGs_A","Carrier_MAGs_B","Reconstruction_error_pp"])
        mdw = writer(stack,output/"Focal_MAG_matched_contributions.tsv.gz",IDFIELDS+["GRESP","Capability","MAG_ID","Domain","Phylum","Order",
                     "Completeness","Contamination","Detected_A","Detected_B","Contribution_A_pp","Contribution_B_pp","B_minus_A_pp","Carrier_detection_status"])
        for filt, _ in FILTERS:
            for subset in ("All373","HQ90_C5"):
                for typ,stratum,site,label,sa,sb in comparisons:
                    wa,wb = weights[(filt,subset,sa)],weights[(filt,subset,sb)]
                    ids = dict(Filter=filt,Subset=subset,Comparison_type=typ,Stratum=stratum,Site=site,Contrast=label,
                               Sample_A=sa,Sample_B=sb,Treatment_A=metadata[sa]["Treatment"],Treatment_B=metadata[sb]["Treatment"],
                               Habitat_A=metadata[sa]["Habitat"],Habitat_B=metadata[sb]["Habitat"])
                    if not wa or not wb:
                        skipped.append(dict(ids,Reason="No_retained_MAGs_in_one_or_both_samples")); continue
                    for rank in members:
                        # All 37 phylum-level traits for condition comparisons;
                        # habitat contrasts and all order-level drill-downs use fixed focal list.
                        fset = features if rank=="Phylum" and typ=="Condition" else FOCAL
                        for c in fset:
                            components, abscomponents = [0.,0.,0.],[0.,0.,0.]
                            individual = []
                            for t in members[rank]:
                                pa,fa,ya,nta,nca = profile[(filt,subset,sa,rank,t,c)]
                                pb,fb,yb,ntb,ncb = profile[(filt,subset,sb,rank,t,c)]
                                r,m,e,status,error = decompose(pa,fa,pb,fb)
                                max_error = max(max_error,abs(error))
                                for i,v in enumerate((r,m,e)):
                                    components[i] += v; abscomponents[i] += abs(v)
                                row = dict(ids,**taxinfo[(rank,t)],GRESP=c[0],Capability=c[1],Taxon_abundance_A_pct=pa,Taxon_abundance_B_pct=pb,
                                           Within_taxon_carrier_A_fraction=fa if fa is not None else "",Within_taxon_carrier_B_fraction=fb if fb is not None else "",
                                           Taxon_contribution_A_pp=ya,Taxon_contribution_B_pp=yb,Net_taxon_difference_pp=yb-ya,
                                           Representation_pp=r,Within_lineage_mixture_pp=m,One_sided_lineage_pp=e,Taxon_detection_status=status,
                                           Detected_MAGs_A=nta,Detected_MAGs_B=ntb,Carrier_MAGs_A=nca,Carrier_MAGs_B=ncb,Reconstruction_error_pp=error)
                                tw.writerow(row)
                                if c in FOCAL:
                                    gkey = tuple(row[k] for k in ["Filter","Subset","Comparison_type","Stratum","Contrast"]+TAXFIELDS+["GRESP","Capability"])
                                    focal_taxon_groups[gkey].append(row)
                                    if ya>0 or yb>0: individual.append(row)
                            qa,qb = prevalence[(filt,subset,sa,c)],prevalence[(filt,subset,sb,c)]
                            close(sum(components),qb-qa,"taxon components sum to prevalence difference")
                            denom = sum(abscomponents)
                            shared_gross = abscomponents[0]+abscomponents[1]
                            row = dict(ids,Rank=rank,GRESP=c[0],Capability=c[1],Prevalence_A_pct=qa,Prevalence_B_pct=qb,Net_difference_pp=qb-qa,
                                       **dict(zip(COMPONENTS,components)),
                                       Sum_abs_representation_pp=abscomponents[0],Sum_abs_mixture_pp=abscomponents[1],Sum_abs_one_sided_pp=abscomponents[2],
                                       Representation_share_absolute_components_pct=100*abscomponents[0]/denom if denom else "",
                                       Mixture_share_absolute_components_pct=100*abscomponents[1]/denom if denom else "",
                                       One_sided_share_absolute_components_pct=100*abscomponents[2]/denom if denom else "",
                                       Shared_lineage_representation_fraction=abscomponents[0]/shared_gross if shared_gross else "",
                                       Sum_abs_components_pp=denom,Reconstruction_error_pp=sum(components)-(qb-qa))
                            site_rows.append(row)
                            # Top entries by mean carrier contribution, not by a favourable test result.
                            individual.sort(key=lambda z: (-(z["Taxon_contribution_A_pp"]+z["Taxon_contribution_B_pp"])/2,z["Taxon_ID"]))
                            top_rows.extend(individual[:5])
                    for c in FOCAL:
                        changes = []
                        for mag in carriers[c]:
                            if mag not in wa and mag not in wb: continue
                            ya,yb = 100*wa.get(mag,0),100*wb.get(mag,0)
                            changes.append(yb-ya); t = taxonomy[mag]
                            status = "Shared_carrier" if mag in wa and mag in wb else "Detected_B_only" if mag in wb else "Detected_A_only"
                            mdw.writerow(dict(ids,GRESP=c[0],Capability=c[1],MAG_ID=mag,
                                **{k:t[k] for k in ("Domain","Phylum","Order")},Completeness=quality[mag][0],Contamination=quality[mag][1],
                                Detected_A=int(mag in wa),Detected_B=int(mag in wb),Contribution_A_pp=ya,Contribution_B_pp=yb,
                                B_minus_A_pp=yb-ya,Carrier_detection_status=status))
                        close(sum(changes),prevalence[(filt,subset,sb,c)]-prevalence[(filt,subset,sa,c)],"individual MAG reconstruction")
        # Equal-site means are computed from per-sample/per-site products, never mean(P)*mean(F).
        taxon_means = []
        tgfields = ["Filter","Subset","Comparison_type","Stratum","Contrast"]+TAXFIELDS+["GRESP","Capability"]
        for gkey, rr in sorted(focal_taxon_groups.items()):
            if len({r["Site"] for r in rr}) != len(rr): raise ValueError("Duplicate site in taxon summary")
            rec = dict(zip(tgfields,gkey),N_sites=len(rr),Sites=",".join(sorted(r["Site"] for r in rr)),
                **{"Mean_"+k:mean(r[k] for r in rr) for k in ["Taxon_abundance_A_pct","Taxon_abundance_B_pct",
                    "Taxon_contribution_A_pp","Taxon_contribution_B_pp","Net_taxon_difference_pp"]+COMPONENTS},
                N_sites_taxon_detected_both=sum(r["Taxon_detection_status"]=="Detected_both" for r in rr),
                N_sites_net_positive=sum(r["Net_taxon_difference_pp"]>1e-10 for r in rr),
                Mean_abs_taxon_change_pp=mean(abs(r["Net_taxon_difference_pp"]) for r in rr))
            taxon_means.append(rec)
        writer(stack,output/"Focal_taxon_matched_site_means.tsv.gz",list(taxon_means[0])).writerows(taxon_means)
        writer(stack,output/"Matched_site_component_summaries.tsv",list(site_rows[0])).writerows(site_rows)
        if top_rows: writer(stack,output/"Focal_top_taxon_contributors.tsv",tw.fieldnames).writerows(top_rows)
        if skipped: writer(stack,output/"Unevaluable_comparisons.tsv",list(skipped[0])).writerows(skipped)
        groups = defaultdict(list)
        groupfields = ["Filter","Subset","Comparison_type","Stratum","Contrast","Rank","GRESP","Capability"]
        for r in site_rows: groups[tuple(r[k] for k in groupfields)].append(r)
        summary = []
        for key,rows in sorted(groups.items()):
            if len({r["Site"] for r in rows}) != len(rows): raise ValueError("Duplicate site in summary")
            means = {"Mean_"+k:mean(r[k] for r in rows) for k in ["Net_difference_pp"]+COMPONENTS+
                     ["Sum_abs_representation_pp","Sum_abs_mixture_pp","Sum_abs_one_sided_pp"]}
            summary.append(dict(zip(groupfields,key),N_sites=len(rows),Sites=",".join(r["Site"] for r in rows),
                                **means,N_sites_net_positive=sum(r["Net_difference_pp"]>1e-10 for r in rows),
                                N_sites_representation_larger_than_mixture=sum(r["Sum_abs_representation_pp"]>r["Sum_abs_mixture_pp"]+1e-10 for r in rows),
                                N_sites_mixture_larger_than_representation=sum(r["Sum_abs_mixture_pp"]>r["Sum_abs_representation_pp"]+1e-10 for r in rows)))
        writer(stack,output/"Matched_contrast_component_means.tsv",list(summary[0])).writerows(summary)
        # Compare settings only on sites evaluable in ALL four settings for the given contrast.
        sg = defaultdict(dict)
        for r in site_rows:
            if (r["GRESP"],r["Capability"]) not in FOCAL: continue
            k = tuple(r[f] for f in ["Comparison_type","Stratum","Contrast","Rank","GRESP","Capability"])
            sg[k].setdefault((r["Filter"],r["Subset"]),{})[r["Site"]] = r
        sensitivity = []
        settings = [(f,s) for f,_ in FILTERS for s in ("All373","HQ90_C5")]
        for key,setting in sorted(sg.items()):
            common = set.intersection(*(set(setting.get(k,{})) for k in settings))
            if not common: continue
            for f,sub in settings:
                rr = [setting[(f,sub)][s] for s in sorted(common)]
                sensitivity.append(dict(zip(["Comparison_type","Stratum","Contrast","Rank","GRESP","Capability"],key),
                    Filter=f,Subset=sub,N_matched_sites=len(common),Sites=",".join(sorted(common)),
                    **{"Mean_"+k:mean(r[k] for r in rr) for k in ["Net_difference_pp"]+COMPONENTS+
                       ["Sum_abs_representation_pp","Sum_abs_mixture_pp","Sum_abs_one_sided_pp"]}))
        if sensitivity: writer(stack,output/"Focal_sensitivity_same_sites.tsv",list(sensitivity[0])).writerows(sensitivity)
    manifest = dict(catalogue_domains=dict(domains),n_MAGs=len(mags),n_metagenomes=len(samples),
                    strict_quality_n=len(hq),strict_quality_domains=dict(Counter(taxonomy[m]["Domain"] for m in hq)),
                    primary_filter="Breadth10",focal_capabilities=FOCAL,checks=checks,max_product_reconstruction_error_pp=max_error,
                    python=sys.version,inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(sources))})
    (output/"Run_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (output/"README.txt").write_text(README)
    print("\n=== VALIDATION ===")
    for s in checks: print(s)
    print("Maximum lineage product reconstruction error:",max_error)
    print("\n=== PHYLUM ACCOUNTING: BESE MINUS BARE, RZ ===")
    print("Signed components in percentage points; R=lineage representation; M=within-lineage carrier mixture; E=one-sided detection.")
    for r in summary:
        if r["Rank"]=="Phylum" and r["Contrast"]=="BESE_minus_BARE" and r["Stratum"]=="RZ" and (r["GRESP"],r["Capability"]) in FOCAL:
            print("{Filter} {Subset} {Capability}: sites={N_sites} net={Mean_Net_difference_pp:+.2f}; R={Mean_Representation_pp:+.2f}; M={Mean_Within_lineage_mixture_pp:+.2f}; E={Mean_One_sided_lineage_pp:+.2f}".format(**r))
    print("\n=== MAIN PHYLUM CONTRIBUTORS TO FOCAL BESE-BARE RZ DIFFERENCES ===")
    print("All373 / Breadth10; ranked by mean absolute site-level contribution difference (not p-values).")
    for c in FOCAL:
        rr = [r for r in taxon_means if r["Filter"]=="Breadth10" and r["Subset"]=="All373" and r["Rank"]=="Phylum"
              and r["Contrast"]=="BESE_minus_BARE" and r["Stratum"]=="RZ" and (r["GRESP"],r["Capability"])==c]
        rr.sort(key=lambda r: (-r["Mean_abs_taxon_change_pp"],r["Taxon_ID"]))
        print(c[1])
        for r in rr[:3]:
            print("  {Domain} / {Phylum}: net={Mean_Net_taxon_difference_pp:+.2f}; R={Mean_Representation_pp:+.2f}; M={Mean_Within_lineage_mixture_pp:+.2f}; E={Mean_One_sided_lineage_pp:+.2f}; sites={N_sites}".format(**r))
    print("\nSUCCESS: Taxonomy-capability abundance decomposition")
    print("All 37 capabilities resolved by phylum; seven focal capabilities resolved to orders and carrier MAGs.")
    print("No p-values calculated: site-matched exact accounting, not an additional null-model screen.")
    print("Output:",output)

README = """LINEAGE-CAPABILITY CONTRIBUTION ACCOUNTING

P_ts = taxon's summed relative mean-depth abundance within all detected MAGs in the selected quality subset.
F_tks = capability-carrying abundance / total abundance of that taxon. F is abundance-weighted, NOT the proportion of genomes by count.
Y_tks = P_ts * F_tks. Output P is percent, F is a fraction, so Y is percentage points of community capability prevalence.
Sum Y across phyla = community prevalence for that capability. Zero-carrier taxa remain in the denominator.

For samples A and B from the same site:
Delta Y = (P_B-P_A)*(F_A+F_B)/2 + (F_B-F_A)*(P_A+P_B)/2.
These are a symmetric arithmetic allocation of change to lineage representation and within-lineage carrier mixture. They do not establish causation.
If the lineage is undetected on one side, F there is undefined. Its entire observed contribution change is retained in a separate one-sided-detection term; no F value is invented. Both-absent lineages contribute zero. No retained MAGs in an entire quality-filtered sample means the sample is unevaluable, NOT zero functional potential.

Condition contrasts are explicitly B minus A: BESE_minus_BARE means BESE - BARE. Habitat contrasts are RZ_minus_BS within the SAME site and condition. Within each habitat, all available site-matched condition comparisons include Natural and Degraded as applicable. No missing metagenomes are imputed. Group means weight available sites equally. Habitat comparisons and sensitivity comparisons retain their matched site sets.

All37 phylum rows and fixed focal order rows retain full GTDB paths, preventing accidental merging of unnamed orders across different phyla/classes. All21 archaeal and all352 bacterial MAGs are included when detected. Individual carrier MAGs have fixed binary annotation calls; within-lineage changes do not imply genes changing inside the same MAG. Broad clades may contain several distinct ecological strategies.

FOCAL SET (chosen from previously discussed contrasts, not these outputs): trehalose; polyphosphate turnover; ammonium uptake; iron-complex uptake; low-oxygen respiration; nitrate respiration; N2O reduction. All other capabilities, including candidate siderophore biosynthesis and both separately keyed DNRA rules, remain in the 37-capability phylum analysis.

The strict quality subset uses estimated completeness >=90% and contamination <5%. Abundances are renormalised within that subset. Sample_catalogue_coverage.tsv reports how much of full-catalogue detected depth remains in each sample at the SAME breadth filter. This is not the fraction of all metagenomic reads recruited. QC is an inclusion sensitivity, not a correction of binary gene absence. No annotation rule is redefined, and archaeal pathway annotation caveats remain.

Sums of absolute component magnitudes quantify contributions before cancellation. Their proportions are NOT fractions of net change, percentages of explained variance or statistical effect attribution. Signed components add to the net prevalence contrast; absolute magnitudes generally do not. Conclusions depend on the taxonomic rank. An order appearing within a phylum present on both sides can be a one-sided order effect but a within-phylum mixture effect. Do not add phylum and order decompositions together: they are alternative resolutions of the same total.

No formal significance tests or new allocation nulls are performed. Individual MAGs, taxa and capabilities are not independent field replicates. Use site effects, both detection filters, quality sensitivities and the original restricted tests when writing inference. Never call a one-sided detection an ecological gain/loss or a specialist without additional evidence.

OUTPUTS
Sample_catalogue_coverage.tsv: denominator, detected richness, domain inclusion and strict-subset retention.
Phylum_capability_profiles.tsv.gz: all samples x phyla x 37 capabilities; P, F, Y and contribution as share of capability.
Focal_order_capability_profiles.tsv.gz: same quantities for every order for the fixed focal set.
Focal_MAG_sample_contributions.tsv.gz: candidate carrier MAGs and sample contributions, mean depth, genome breadth and quality.
Matched_taxon_decomposition.tsv.gz: exact lineage components, status, P and F on each side.
Matched_site_component_summaries.tsv: aggregate effects separately at each rank, one row per capability and site comparison.
Matched_contrast_component_means.tsv: equal-site means, directional consistency and absolute magnitudes (not significance tests).
Focal_top_taxon_contributors.tsv: top five phyla/orders per focal pair by mean carrier contribution; display aid, not additional testing.
Focal_taxon_matched_site_means.tsv.gz: all focal phylum/order equal-site means, including zero contributions; average sample products, not products of averages.
Focal_MAG_matched_contributions.tsv.gz: identities and per-MAG changes, including shared/one-sided carrier labels.
Focal_sensitivity_same_sites.tsv: four breadth/quality settings compared on identical evaluable sites.
Run_manifest.json: input hashes, checks and subset counts.

Compressed .tsv.gz tables can be read directly in R with read.delim(gzfile('filename.tsv.gz')).
"""

if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, FileNotFoundError, OSError) as exc:
        print("ERROR:",exc,file=sys.stderr)
        sys.exit(1)
