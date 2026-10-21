#!/usr/bin/env python3
"""Read-only source-grounded HQ sensitivity + candidate iron evidence audit. Python 3.6 compatible.
Run from project root: python3 scripts/validate_quality_iron.py
"""
import csv
from pathlib import Path
from collections import defaultdict
ROOT=Path.cwd(); OUT=ROOT/'results/ecological_strategies/final_validation';OUT.mkdir(parents=True,exist_ok=True)
def load(path,key):
    with path.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==len(set(x[key] for x in rows))
    return {r[key]:r for r in rows}
def save(name,rows):
    if not rows:return
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
Q=load(ROOT/'results/final_summaries/combined_binners/final_373_MAGs_CheckM2_quality.tsv','Final_MAG')
T=load(ROOT/'results/final_summaries/combined_binners/final_373_MAGs_GTDB_taxonomy.tsv','Final_MAG')
D=load(ROOT/'mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_mean_depth_matrix.tsv','MAG_ID')
B=load(ROOT/'mapping/competitive_MAGs/MAG_abundance/final_373_MAGs_breadth_matrix.tsv','MAG_ID')
assert len(Q)==len(T)==len(D)==len(B)==373 and set(Q)==set(T)==set(D)==set(B)
HQ={m for m,r in Q.items() if float(r['Completeness'])>=90 and float(r['Contamination'])<5}
S=list(next(iter(D.values())));S.remove('MAG_ID');assert len(S)==20
from collections import OrderedDict
spec=OrderedDict([
 ('Osmotic',('osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv','Capability_Trehalose Capability_Compatible_solute_uptake Capability_Ectoine_biosynthesis Capability_Mechanosensitive_protection Capability_Ion_homeostasis')),
 ('Redox',('redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv','Capability_Aerobic_respiration Capability_LowO2_respiration Capability_Nitrate_respiration Capability_Nitrite_to_NO Capability_DNRA Capability_NO_reduction Capability_N2O_reduction Capability_Dissimilatory_sulfate_reduction_core Capability_Reduced_sulfur_oxidation Capability_Sulfite_oxidation')),
 ('Nitrogen',('nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv','Capability_Nitrogen_fixation Capability_Ammonium_uptake Capability_Assimilatory_nitrate_reduction Capability_Urea_utilisation Capability_DNRA Capability_Complete_denitrification Capability_Nitrification Capability_Anammox')),
 ('Carbon_P',('carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv','Capability_Carbon_fixation Capability_C1_utilisation Capability_High_affinity_phosphate_uptake Capability_Organic_P_mobilisation Capability_Phosphonate_acquisition Capability_PolyP_turnover')),
 ('Iron',('iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv','Ferrous_iron_uptake Direct_ferric_iron_uptake Siderophore_biosynthesis Iron_complex_uptake')),
 ('Colonisation',('colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv','Chemotaxis Motility Attachment_biofilm_matrix Direct_plant_interaction'))])
cap={}; calls=[]
for group,(rel,cols) in spec.items():
  a=load(ROOT/'results/ecological_strategies'/rel,'Final_MAG');assert set(a)==set(Q)
  for col in cols.split():
    carrier={m for m,r in a.items() if float(r[col])==1}
    assert all(float(r[col]) in (0,1) for r in a.values())
    cap[(group,col)]=carrier
    calls.append({'GRESP':group,'Capability':col,'All_carriers':len(carrier),'All_MAGs':373,'HQ_carriers':len(carrier & HQ),'HQ_MAGs':len(HQ),'All_pct':round(100*len(carrier)/373,4),'HQ_pct':round(100*len(carrier & HQ)/len(HQ),4)})
assert len(cap)==37
save('HQ_37capability_catalogue.tsv',calls)
rows=[]
for s in S:
  parts=s.split('_');site,tr,hab=parts[1:4]
  for filt,cut in (('Breadth10',.1),('Breadth50',.5)):
    retain={m:float(D[m][s]) for m in Q if float(D[m][s])>0 and float(B[m][s])>=cut}
    for label,subset in (('All373',set(Q)),('HQ90_C5',HQ)):
      w={m:v for m,v in retain.items() if m in subset}; tot=sum(w.values())
      for (g,c),carriers in cap.items():
        rows.append({'Sample':s,'Site':site,'Treatment':tr,'Habitat':hab,'Filter':filt,'Subset':label,'GRESP':g,'Capability':c,'Detected_MAGs':len(w),'Retained_total_depth':round(tot,8),'Prevalence_pct':round(100*sum(v for m,v in w.items() if m in carriers)/tot,6) if tot>0 else ''})
save('HQ_37capability_sample_prevalence.tsv',rows)
producers=cap[('Iron','Siderophore_biosynthesis')]
prod=[]
for m in sorted(producers):
  tax=T[m]['GTDB_classification'];g=next((v[3:] for v in tax.split(';') if v.startswith('g__')),'')
  for s in S:
    prod.append({'MAG_ID':m,'GTDB_genus':g,'Completeness':Q[m]['Completeness'],'Contamination':Q[m]['Contamination'],'HQ90_C5':int(m in HQ),'Sample':s,'Site':s.split('_')[1],'Treatment':s.split('_')[2],'Habitat':s.split('_')[3],'Mean_depth':D[m][s],'Breadth':B[m][s],'Detected_Breadth10':int(float(D[m][s])>0 and float(B[m][s])>=.1),'Detected_Breadth50':int(float(D[m][s])>0 and float(B[m][s])>=.5)})
save('Iron_three_producers_20samples.tsv',prod)
master=ROOT/'results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv'
if master.exists():
  from collections import defaultdict
  lookup=defaultdict(list)
  with master.open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
      if row.get('Final_MAG') in producers:lookup[row['Final_MAG']].append(row)
  routes={'Aerobactin':['K03897','K03896','K03894','K03895'],'Mycobactin':['K04781','K04787','K04788','K04789','K04790','K04791','K04792','K04793']}
  evid=[]
  for m in sorted(producers):
    for route,kos in routes.items():
      for ko in kos:
        found=[r for r in lookup[m] if r.get('KO')==ko]
        evid.append({'MAG_ID':m,'Route':route,'KO':ko,'Detected_genes':len(found),'Gene_IDs':';'.join(x.get('gene_id','') for x in found),'KO_definitions':';'.join(sorted(set(x.get('KO_definition','') for x in found)))})
  save('Iron_aerobactin_mycobactin_gene_evidence.tsv',evid)
  print('Gene evidence available:',sum(len(v) for v in lookup.values()),'KO hits among 3 candidates')
else:print('GENE MASTER NOT FOUND: locus-level verification requires Gadi gene-master file')
print('SUCCESS: HQ and iron candidate audit')
print('HQ MAGs: {}/373'.format(len(HQ)))
print('Bacterial HQ:',sum(T[m]['GTDB_domain']=='Bacteria' for m in HQ))
print('Siderophore candidate MAGs:',len(producers),'HQ:',len(producers & HQ))
print('Outputs:',OUT)
