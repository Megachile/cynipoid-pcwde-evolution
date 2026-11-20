"""Author-approved Table S3 protein-count substitutions, with an audit per numeral."""
import re

def project(text, numbers, document):
    audit = []
    def put(pattern, key):
        nonlocal text
        value = numbers
        for part in key.split('.'):
            value = value[part]
        if not isinstance(value, int):
            raise RuntimeError(f'outgroup protein count is not integer: {key}')
        def sub(match):
            audit.append(dict(document=document, key=key, old=match.group('count'), new=str(value)))
            start = match.start('count') - match.start()
            end = match.end('count') - match.start()
            return match.group()[:start] + str(value) + match.group()[end:]
        text, count = re.subn(pattern, sub, text)
        if count != 1:
            raise RuntimeError(f'outgroup count source drift for {key}: {count}')
    if document == 'manuscript':
        put(r'(?P<count>\d+) unique PCWDE outgroup sequences \(Table S3\)', 'supplement.outgroup_parent_proteins_total')
    else:
        put(r'(?P<count>\d+) total PCWDE outgroup sequences \(Table S3\)', 'supplement.outgroup_parent_proteins_total')
        for family in ('PL1_8', 'PL4_4', 'GH5_2'):
            escaped = re.escape(family.replace('_', '\\_'))
            prefix = r'For ' + escaped + r', we found '
            put(prefix + r'(?P<count>\d+) bacterial sequences', 'supplement.' + family + '.screen_bacterial_parent_proteins')
            put(prefix + r'\d+ bacterial sequences and (?P<count>\d+) non-cynipoid eukaryote sequences', 'supplement.' + family + '.screen_eukaryotic_parent_proteins')
        for family, taxa in (('PL1_8', ('midge','bee','oomycete')), ('PL4_4', ('fungus','plant','beetle')), ('GH5_2', ('nematode','protist','beetle'))):
            # Restrict each replacement to its own sentence to disambiguate taxa.
            start = text.index('For ' + family.replace('_', '\\_') + ', we found ')
            end = text.index(').', start) + 2
            sentence = text[start:end]
            labels = {'midge':'Cecidomyiidae','bee':'Apidae','oomycete':'Oomycota','fungus':'fungus','plant':r'\*Cuscuta\*','beetle':'Coleoptera','nematode':'Nematoda','protist':'Protista'}
            for taxon in taxa:
                key='supplement.'+family+'.outgroup_parent_taxa.'+taxon
                value=numbers['supplement'][family]['outgroup_parent_taxa'][taxon]
                pattern=r'(?P<count>\d+) '+labels[taxon]
                hits=list(re.finditer(pattern,sentence))
                if len(hits)!=1:raise RuntimeError(f'outgroup taxon source drift: {key}')
                match=hits[0];audit.append(dict(document=document,key=key,old=match.group('count'),new=str(value)))
                sentence=sentence[:match.start('count')]+str(value)+sentence[match.end('count'):]
            text=text[:start]+sentence+text[end:]
    return text, audit
