"""Regress verified transient cleanup, source mutation and confined paths."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lifecycle import inventory, consume_source, consume_comparison

with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp); repo=root/'repo';repo.mkdir();receipt=root/'receipt';receipt.mkdir()
    source=root/'staged_input_fixture';source.mkdir();(source/'payload').write_bytes(b'original\n')
    stage=root/'copied';shutil.copytree(source,stage)
    expected=inventory(source);seal=root/'inventory.json';seal.write_text(json.dumps(expected))
    spec={'source':str(source),'inventory':str(seal),'inventory_sha256':hashlib.sha256(seal.read_bytes()).hexdigest()}
    cfg={'bounded_staging':spec,'destination':str(root/'output'),'receipt':str(receipt)}
    cases=[]
    def reject(name,call,token):
        try:call()
        except ValueError as e:assert str(e)==token,(name,str(e));cases.append(name)
        else:raise AssertionError(name+' accepted')
        assert source.exists()
    (stage/'payload').write_bytes(b'altered\n')
    reject('copy mutation',lambda:consume_source(cfg,repo,source,stage,expected,receipt),'SOURCE_CLEANUP_IDENTITY_FAIL')
    (stage/'payload').write_bytes(b'original\n')
    (source/'payload').write_bytes(b'altered\n')
    reject('source mutation',lambda:consume_source(cfg,repo,source,stage,expected,receipt),'SOURCE_CLEANUP_IDENTITY_FAIL')
    (source/'payload').write_bytes(b'original\n')
    spec['inventory_sha256']='0'*64
    reject('seal mutation',lambda:consume_source(cfg,repo,source,stage,expected,receipt),'SOURCE_CLEANUP_SEAL_FAIL')
    spec['inventory_sha256']=hashlib.sha256(seal.read_bytes()).hexdigest()
    spec['source']=str(root/'good_data')
    reject('wrong source',lambda:consume_source(cfg,repo,source,stage,expected,receipt),'SOURCE_CLEANUP_CONFINEMENT_FAIL')
    spec['source']=str(source)
    consume_source({},repo,source,stage,expected,receipt);assert source.exists();cases.append('default preserves source')
    consume_source(cfg,repo,source,stage,expected,receipt);assert not source.exists() and inventory(stage)==expected;cases.append('verified source consumption')
    comparison=root/'output.stage';shutil.copytree(stage,comparison);log=root/'stage.log';log.write_text('STAGE_COMPARE_FAIL files=1 differing=1\n')
    try:consume_comparison(comparison,repo,inventory(comparison),log,receipt)
    except ValueError as e:assert str(e)=='COMPARE_CLEANUP_GATE_FAIL'
    else:raise AssertionError('failed compare deleted')
    assert comparison.exists();cases.append('failed compare preserved')
    log.write_text('STAGE_COMPARE_PASS files=1 differing=0\n')
    consume_comparison(comparison,repo,inventory(comparison),log,receipt);assert not comparison.exists();cases.append('validated comparison consumption')
    print(json.dumps({'status':'PASS','cases':cases}))
