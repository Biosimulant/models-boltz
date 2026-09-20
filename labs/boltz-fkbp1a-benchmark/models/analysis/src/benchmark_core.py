"""Shared scientific transformations; no platform integration or hidden predictions."""
from __future__ import annotations
import hashlib,io,json,math,statistics
from pathlib import Path

def descriptors(bundle):
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    result=[]
    for ligand in bundle['ligands']:
        molecule=Chem.MolFromSmiles(ligand['smiles'])
        if molecule is None:raise ValueError('Invalid isomeric SMILES: '+ligand['name'])
        result.append({**ligand,'canonical_isomeric_smiles':Chem.MolToSmiles(molecule,isomericSmiles=True),
                       'molecular_mass_g_mol':Descriptors.MolWt(molecule),
                       'heavy_atom_count':molecule.GetNumHeavyAtoms()})
    return result

def check_affinity(affinity):
    score=affinity['affinity_pred_value'];prob=affinity['affinity_probability_binary']
    for value in (score,prob):
        if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value):
            raise ValueError('Non-finite or non-numeric scientific output')
    if not 0<=prob<=1:raise ValueError('Binding probability outside [0,1]')
    return {'affinity_pred_value':score,'affinity_probability_binary':prob,
            'derived_p_scale':6-score,'affinity_interpretation':'predicted mixed-label affinity/activity score; not measured IC50 or Kd'}

def aligned_ca_rmsd(predicted,reference):
    import numpy as np
    from Bio import Align
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    from Bio.SVDSuperimposer import SVDSuperimposer
    def chains(text):
        structure=MMCIFParser(QUIET=True,auth_chains=False).get_structure('structure',io.StringIO(text))
        return [(chain.id,[r for r in chain if r.id[0]==' ' and 'CA' in r]) for chain in structure[0]]
    pred=max(chains(predicted),key=lambda x:len(x[1]))
    refs=[c for c in chains(reference) if len(c[1])>0]
    aligner=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-5,extend_gap_score=-.5)
    ps=''.join(seq1(r.resname) for r in pred[1])
    best=None
    for ref in refs:
        rs=''.join(seq1(r.resname) for r in ref[1]);alignment=aligner.align(ps,rs)[0]
        if best is None or alignment.score>best[0]:best=(alignment.score,ref,rs,alignment)
    if best is None:raise ValueError('No reference protein C-alpha residues')
    _,ref,rs,alignment=best
    pc=[];rc=[]
    for pblock,rblock in zip(*alignment.aligned):
        for i,j in zip(range(*pblock),range(*rblock)):
            if ps[i]==rs[j] and ps[i]!='X':
                pc.append(pred[1][i]['CA'].coord);rc.append(ref[1][j]['CA'].coord)
    if len(pc)<3:raise ValueError('Insufficient sequence-matched C-alpha atoms')
    sup=SVDSuperimposer();sup.set(np.asarray(rc),np.asarray(pc));sup.run()
    return {'ca_rmsd_angstrom':float(sup.get_rms()),'aligned_exact_residues':len(pc),
            'predicted_ca_count':len(pred[1]),'reference_ca_count':len(ref[1]),
            'predicted_alignment_coverage':len(pc)/len(pred[1]),'reference_chain':ref[0],
            'reference':'PDB 1FKJ','reference_sha256':hashlib.sha256(reference.encode()).hexdigest()}

def geometry(cif,expected_heavy_atoms):
    import numpy as np
    from Bio.PDB import MMCIFParser
    structure=MMCIFParser(QUIET=True,auth_chains=False).get_structure('prediction',io.StringIO(cif))[0]
    protein=[a for a in structure['A'].get_atoms() if a.element not in ('H','D')]
    ligand=[a for a in structure['B'].get_atoms() if a.element not in ('H','D')]
    if len(ligand)!=expected_heavy_atoms:raise ValueError(f'Ligand atom count mismatch: {len(ligand)} != {expected_heavy_atoms}')
    if not protein or not ligand:raise ValueError('Missing protein or ligand atoms')
    p=np.array([a.coord for a in protein]);l=np.array([a.coord for a in ligand])
    if not np.isfinite(p).all() or not np.isfinite(l).all():raise ValueError('Non-finite coordinates')
    d=np.linalg.norm(p[:,None,:]-l[None,:,:],axis=-1)
    close=d<=4.0
    residues={protein[i].get_parent().id[1] for i in np.where(close.any(axis=1))[0]}
    return {'protein_heavy_atoms':len(protein),'ligand_heavy_atoms':len(ligand),
            'contact_atom_pairs_4A':int(close.sum()),'contact_residue_count_4A':len(residues),
            'contact_residue_indices_4A':sorted(residues),'minimum_heavy_atom_distance_angstrom':float(d.min()),
            'contact_cutoff_angstrom':4.0,'contact_interpretation':'geometric proximity only; not hydrogen-bond assignment'}

def summarize(rows,question):
    expected={1:1,2:3,3:3}[question]
    completed=[r for r in rows if r.get('status')=='completed']
    result={'question':question,'expected_count':expected,'completed_count':len(completed),'rows':rows,
            'status':'completed' if len(completed)==expected and len(rows)==expected else 'partial',
            'units':{'molecular_mass_g_mol':'g/mol','ca_rmsd_angstrom':'angstrom','affinity_pred_value':'log10(micromolar concentration scale); mixed biochemical labels','affinity_probability_binary':'1'}}
    if question==2:
        result['predicted_score_ranking']=[r['ligand'] for r in sorted(completed,key=lambda r:(r['affinity_pred_value'],r['ligand']))]
        result['ranking_claim']='Exploratory model score ordering; no experimental best-ligand claim'
    if question==3:
        result['statistics']={}
        for key in ('affinity_pred_value','affinity_probability_binary','ca_rmsd_angstrom'):
            values=[r[key] for r in completed if key in r]
            result['statistics'][key]={'n':len(values),'mean':statistics.mean(values) if values else None,
                'sample_sd':statistics.stdev(values) if len(values)>1 else None,
                'min':min(values) if values else None,'max':max(values) if values else None}
    return result

def analyze_native(root,bundle,jobs,reference):
    root=Path(root);ds={d['name']:d for d in descriptors(bundle)};rows=[]
    for job in jobs:
        base=root/job['id'];row={**job,**ds[job['ligand']]}
        try:
            affinities=list(base.rglob('affinity_*.json'));structures=list(base.rglob('*model_0.cif'));confidence=list(base.rglob('confidence_*model_0.json'))
            if len(affinities)!=1 or len(structures)!=1 or len(confidence)!=1:raise ValueError('Missing or ambiguous primary artifacts')
            af=json.loads(affinities[0].read_text());cf=json.loads(confidence[0].read_text());cif=structures[0].read_text()
            row.update(check_affinity(af));row.update(geometry(cif,row['heavy_atom_count']))
            if job['question'] in (1,3):row.update(aligned_ca_rmsd(cif,reference))
            row.update(status='completed',confidence=cf,structure_sha256=hashlib.sha256(cif.encode()).hexdigest())
        except Exception as error:row.update(status='failed',error=str(error))
        rows.append(row)
    return {str(q):summarize([r for r in rows if r['question']==q],q) for q in sorted({j['question'] for j in jobs})}
