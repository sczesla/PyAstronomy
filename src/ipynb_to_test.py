from pathlib import Path
from nbconvert import ScriptExporter
import nbformat


def insert_mpl_agg(s):
    # Find last import ...
    x = 0
    for i, l in enumerate(s):
        if l.find("import") != -1:
            x = i
    # ... and import matplotlib to set backend to Agg
    s.insert(x+1,"matplotlib.use('Agg')")
    s.insert(x+1,"import matplotlib")
    return s

def del_future(s):
    # Remove __future__ import
    for i, l in enumerate(s):
        if l.find("from __future__") != -1:
            s.pop(i)
            return s
    return s


def sanitize(s, name="example"):
    x = """
import unittest
class TestSanityOfNotebook(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass
        
    def testsanity_%s(self):""" % name
    x = x.split("\n")
    s = [" "*8+l for l in s]
    x.extend(s)
    return x


fns = list(Path(".").rglob("ex_*.ipynb"))

print("Files matching criteria:")
for i, fn in enumerate(fns,1):
    print(i, fn)
print()

for fn in fns:
    if str(fn).find("_build") != -1:
        continue
    c = ''.join(open(str(fn), 'rt').readlines())
    nb = nbformat.reads(c, as_version=4)
    se = ScriptExporter()
    (body, resources) = se.from_notebook_node(nb)
    s = body.split("\n")
    
    s = del_future(s)
    s = insert_mpl_agg(s)
    s = sanitize(s)
    
    ofp = fn.with_name("sanity_"+fn.stem+".py")
    print("Input notebook: ", fn)
    print("Output test: ", ofp)
    
    open(str(ofp), 'wt').writelines([l+"\n" for l in s])
    
