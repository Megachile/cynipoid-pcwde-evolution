"""Small Newick reader and Fitch root-state probe for the fixed species topology."""
from __future__ import annotations
from dataclasses import dataclass, field

@dataclass
class Node:
    name: str = ''
    children: list = field(default_factory=list)
    length: float | None = None

def parse_newick(source):
    tokens=[]; i=0
    while i<len(source):
        if source[i] in '(),:;':tokens.append(source[i]);i+=1
        elif source[i].isspace():i+=1
        else:
            j=i
            while j<len(source) and source[j] not in '(),:; \t\n':j+=1
            tokens.append(source[i:j]);i=j
    pos=0
    def node():
        nonlocal pos
        children=[]
        if tokens[pos]=='(':
            pos+=1
            while True:
                children.append(node())
                if tokens[pos]==',':pos+=1;continue
                if tokens[pos]==')':pos+=1;break
                raise ValueError(('bad Newick',pos,tokens[pos]))
        name=''
        if pos<len(tokens) and tokens[pos] not in ',):;':name=tokens[pos];pos+=1
        length=None
        if pos<len(tokens) and tokens[pos]==':':pos+=1;length=float(tokens[pos]);pos+=1
        return Node(name,children,length)
    result=node()
    if tokens[pos]!=';':raise ValueError('Newick does not end at semicolon')
    return result

def prune(root, keep):
    """Return induced rooted tree; collapse degree-one paths."""
    if not root.children:return Node(root.name,[],root.length) if root.name in keep else None
    children=[x for child in root.children if (x:=prune(child,keep)) is not None]
    if not children:return None
    if len(children)==1:
        child=children[0]
        if root.length is not None and child.length is not None:child.length+=root.length
        return child
    return Node(root.name,children,root.length)

def fitch(root, ones):
    """Return possible states at each node and minimum equal-cost changes."""
    states={};changes=0
    def visit(node):
        nonlocal changes
        if not node.children:
            result={1 if node.name in ones else 0}
        else:
            parts=[visit(child) for child in node.children]
            result=set(parts[0])
            for part in parts[1:]:
                common=result & part
                if common:result=common
                else:result|=part;changes+=1
        states[id(node)]=result
        return result
    visit(root)
    return states,changes
