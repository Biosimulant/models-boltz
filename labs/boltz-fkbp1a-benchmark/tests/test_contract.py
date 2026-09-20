import math
import pytest
from benchmark_core import check_affinity,descriptors,summarize,aligned_ca_rmsd

@pytest.mark.parametrize('x',[float('nan'),float('inf'),True])
def test_invalid_numeric_score_is_rejected(x):
    with pytest.raises(ValueError):check_affinity({'affinity_pred_value':x,'affinity_probability_binary':.5})

def test_log_scale_conversion_and_probability_bounds():
    assert check_affinity({'affinity_pred_value':-3.,'affinity_probability_binary':1.})['derived_p_scale']==9.
    with pytest.raises(ValueError):check_affinity({'affinity_pred_value':0.,'affinity_probability_binary':1.01})

def test_ranking_direction_and_no_imputed_failed_scores():
    rows=[{'ligand':'a','status':'completed','affinity_pred_value':2.},{'ligand':'b','status':'completed','affinity_pred_value':-1.},{'ligand':'c','status':'failed'}]
    result=summarize(rows,2)
    assert result['predicted_score_ranking']==['b','a']
    assert result['status']=='partial'

def test_sample_sd_is_derived_and_missing_is_not_zero():
    r=summarize([{'status':'completed','affinity_pred_value':v,'affinity_probability_binary':.5} for v in (1.,2.,3.)],3)
    assert r['statistics']['affinity_pred_value']['sample_sd']==1.
    assert r['statistics']['ca_rmsd_angstrom']['mean'] is None

def test_stereochemistry_survives_descriptors():
    bundle={'ligands':[{'name':'L','smiles':'N[C@@H](C)C(=O)O'},{'name':'D','smiles':'N[C@H](C)C(=O)O'}]}
    a,b=descriptors(bundle)
    assert a['molecular_mass_g_mol']==b['molecular_mass_g_mol']
    assert a['canonical_isomeric_smiles']!=b['canonical_isomeric_smiles']
