"""Focused counting-unit and multi-genome membership tests."""
from comprehensive_coverage import comprehensive_coverage
R='NODE_38372_length_2218_cov_4.820177:1323-346'
copies=[dict(locus_id=k,family_id='PL1_8',species_id=s) for k,s in [('a','A'),('b','A'),('c','B'),(R,'C')]]
members=[dict(accepted_locus_id=k,family='PL1_8',block_id='one',assignment_status='assigned') for k in ['a','b','c']]
states=[dict(family='PL1_8',block_id=b,genome=g,status='PRESENT',target_unit_ids=ids) for b,g,ids in [('one','A','a;b'),('one','B','c'),('two','A','a;b'),('two','B','c')]]
e,t,c,m=comprehensive_coverage(copies,members,states,{'A':'A','B':'B'})
assert t['PL1_8']==4 and c['PL1_8']==3 and len(m)==2
assert next(r for r in e if r['locus_id']==R)['included']==0
# Two tandem copies in just one genome never establish a multi-genome block.
e,t,c,m=comprehensive_coverage(copies,[r for r in members if r['accepted_locus_id']!='c']+[dict(accepted_locus_id='c',family='PL1_8',block_id='',assignment_status='unassigned')],states[0:1],{'A':'A'})
assert not m and not c
try:
 comprehensive_coverage(copies+copies[:1],members,states,{'A':'A','B':'B'})
except ValueError:
 pass
else:
 raise AssertionError('duplicate unit admitted')
print('COMPREHENSIVE_COUNTING_TEST_PASS')
