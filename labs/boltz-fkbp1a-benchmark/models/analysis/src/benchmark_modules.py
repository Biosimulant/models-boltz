"""Finite BioSimulant composition of chemistry, Boltz inference and structural analysis."""
from __future__ import annotations
import hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
from biosim import BioModule,ExecutionPolicy
from biosim.signals import SignalSpec,AcceptedSignalProfile,unwrap_payload
from .benchmark_core import descriptors,check_affinity,geometry,aligned_ca_rmsd,summarize

SCHEMA={'payload':'json'}
def output_spec(description):
    return SignalSpec.record(schema=SCHEMA,description=description)
def input_spec(description):
    return SignalSpec.record(schema=SCHEMA,accepted_profiles=(AcceptedSignalProfile(signal_type='record',schema=SCHEMA),),description=description)

class CandidateModel(BioModule):
    execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self,bundle=None,jobs=None,question=1):
        self.bundle=bundle;self.jobs=jobs;self.question=int(question)
    def inputs(self):return {}
    def outputs(self):return {'candidates':output_spec('Molecular identities, descriptors with explicit per-field units, sequence, query MSA and run settings.')}
    def execute(self,inputs,*,context):
        if self.question not in (1,2,3):raise ValueError('question must be 1, 2 or 3')
        if not self.bundle or not self.jobs:raise ValueError('Missing frozen input bundle or jobs')
        query=''.join(x for x in self.bundle['msa'].splitlines() if not x.startswith('>'))
        if query!=self.bundle['sequence']:raise ValueError('MSA query mismatch')
        ds=descriptors(self.bundle)
        return {'candidates':{'bundle':self.bundle,'jobs':[j for j in self.jobs if j['question']==self.question],
                'descriptors':ds,'question':self.question,'field_units':{'molecular_mass_g_mol':'g/mol','heavy_atom_count':'1'}}}

class PredictionModel(BioModule):
    execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self,cache_dir=None,command_timeout_s=1800):
        self.cache_dir=cache_dir or os.environ.get('BOLTZ_CACHE','/tmp/boltz-benchmark-cache')
        self.command_timeout_s=int(command_timeout_s)
    def inputs(self):return {'candidates':input_spec('Frozen candidate identities and inference controls.')}
    def outputs(self):return {'predictions':output_spec('Raw structure text, confidence and affinity with field-level physical conventions, plus execution evidence.')}
    def execute(self,inputs,*,context):
        import yaml
        d=unwrap_payload(inputs['candidates']);bundle=d['bundle'];rows=[]
        root=Path(tempfile.mkdtemp(prefix='boltz-benchmark-'))
        (root/'query.a3m').write_text(bundle['msa'])
        for job in d['jobs']:
            ligand=next(l for l in d['descriptors'] if l['name']==job['ligand'])
            jp=root/job['id'];jp.mkdir()
            request={'version':1,'sequences':[{'protein':{'id':'A','sequence':bundle['sequence'],'msa':str(root/'query.a3m')}},
                     {'ligand':{'id':'B','smiles':ligand['smiles']}}],'properties':[{'affinity':{'binder':'B'}}]}
            (jp/'request.yaml').write_text(yaml.safe_dump(request))
            cmd=['boltz','predict',str(jp/'request.yaml'),'--out_dir',str(jp/'outputs'),'--cache',self.cache_dir,
                 '--accelerator','gpu','--devices','1','--output_format','mmcif','--seed',str(job['seed']),
                 '--no_trifast']
            for k,v in bundle['settings'].items():cmd.extend(['--'+k,str(v)])
            row={**job,**ligand,'command':cmd,'attention_backend':'upstream PyTorch fallback (--no_trifast); T4 kernel compatibility'};start=time.monotonic()
            try:
                with (jp/'inference.log').open('w') as log:
                    process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=self.command_timeout_s)
                row['exit_code']=process.returncode
                if process.returncode:raise RuntimeError('Boltz nonzero exit: '+(jp/'inference.log').read_text()[-3000:])
                structures=list(jp.rglob('*model_0.cif'));affinities=list(jp.rglob('affinity_*.json'));confidences=list(jp.rglob('confidence_*model_0.json'))
                if len(structures)!=1 or len(affinities)!=1 or len(confidences)!=1:raise ValueError('Missing or ambiguous primary output artifacts')
                cif=structures[0].read_text();affinity=json.loads(affinities[0].read_text())
                row.update(check_affinity(affinity))
                row.update(status='completed',structure_mmcif=cif,raw_affinity=affinity,confidence=json.loads(confidences[0].read_text()),
                           structure_sha256=hashlib.sha256(cif.encode()).hexdigest())
            except subprocess.TimeoutExpired:row.update(status='timeout',error='Prediction exceeded bounded runtime')
            except Exception as e:row.update(status='failed',error=str(e))
            row['wall_seconds']=time.monotonic()-start;rows.append(row)
            if row['status']!='completed':break
        for job in d['jobs'][len(rows):]:rows.append({**job,'status':'not_run_after_failure'})
        from importlib.metadata import distributions
        versions='\n'.join(sorted(f"{p.metadata['Name']}=={p.version}" for p in distributions() if p.metadata.get('Name')))
        try:
            gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader'],text=True,timeout=15).strip()
        except Exception as error:
            gpu='unavailable: '+type(error).__name__
        checkpoints={}
        for path in Path(self.cache_dir).glob('*.ckpt'):
            h=hashlib.sha256()
            with path.open('rb') as f:
                for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
            checkpoints[path.name]={'sha256':h.hexdigest(),'bytes':path.stat().st_size}
        return {'predictions':{'question':d['question'],'rows':rows,'packages':versions,'gpu':gpu,'checkpoints':checkpoints,
                              'field_units':{'structure_mmcif':'Cartesian angstrom coordinates','affinity_pred_value':'native log10 micromolar scale; mixed labels','affinity_probability_binary':'1'}}}

class AnalysisModel(BioModule):
    execution_policy=ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self,reference_sha256=None):self.reference_sha256=reference_sha256
    def inputs(self):return {'predictions':input_spec('Predictions preserve raw structure/score identities and units.')}
    def outputs(self):return {'report':output_spec('Terminal mixed-unit report: coordinates/RMSD in angstrom, mass in g/mol, probability dimensionless; native affinity convention documented.')}
    def execute(self,inputs,*,context):
        import httpx
        d=unwrap_payload(inputs['predictions']);reference=None
        if d['question'] in (1,3):
            r=httpx.get('https://files.rcsb.org/download/1FKJ.cif',timeout=45);r.raise_for_status()
            if hashlib.sha256(r.content).hexdigest()!=self.reference_sha256:raise ValueError('Experimental reference checksum changed')
            reference=r.text
        rows=[]
        for original in d['rows']:
            row=dict(original)
            if row['status']=='completed':
                try:
                    cif=row['structure_mmcif']
                    row.update(geometry(cif,row['heavy_atom_count']))
                    if reference is not None:row.update(aligned_ca_rmsd(cif,reference))
                except Exception as e:row.update(status='failed_analysis',error=str(e))
            rows.append(row)
        report=summarize(rows,d['question'])
        report.update(gpu=d['gpu'],packages=d['packages'],checkpoints=d['checkpoints'])
        return {'report':report}
