"""Runtime transport regression; simulated CLI failure, never benchmark inference."""
import json,sys
from pathlib import Path
import yaml
from biosim import BioWorld
from biosim.signals import unwrap_payload

LAB=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(LAB/'models/candidates'))
from src.benchmark_modules import CandidateModel,PredictionModel,AnalysisModel

def test_real_world_graph_unwraps_records_and_preserves_failure(monkeypatch,tmp_path):
    import src.benchmark_modules as modules
    parameters={x['name']:x['value'] for x in yaml.safe_load((LAB/'models/candidates/model.yaml').read_text())['parameters']}
    bundle=parameters['bundle'];jobs=parameters['jobs']
    calls=[]
    def simulated_failure(cmd,stdout,**kwargs):
        request=yaml.safe_load(Path(cmd[2]).read_text())
        assert request['sequences'][0]['protein']['sequence']==bundle['sequence']
        assert request['sequences'][1]['ligand']['smiles']==bundle['ligands'][0]['smiles']
        calls.append(cmd)
        stdout.write('intentional test-only CLI failure')
        return type('Exit',(),{'returncode':99})()
    monkeypatch.setattr(modules.subprocess,'run',simulated_failure)
    monkeypatch.setattr(modules.subprocess,'check_output',lambda *a,**kw:'test-only provenance fixture')
    world=BioWorld(communication_step=.01)
    world.add_biomodule('candidates',CandidateModel(bundle,jobs,2))
    world.add_biomodule('prediction',PredictionModel(str(tmp_path)))
    world.add_biomodule('analysis',AnalysisModel())
    world.connect('candidates.candidates','prediction.candidates')
    world.connect('prediction.predictions','analysis.predictions')
    world.run(duration=.01)
    report=unwrap_payload(world.get_outputs('analysis')['report'])
    assert len(calls)==1
    assert report['question']==2 and report['completed_count']==0
    assert report['rows'][0]['exit_code']==99
    assert [r['status'] for r in report['rows']]==['failed','not_run_after_failure','not_run_after_failure']
    assert report['predicted_score_ranking']==[]
