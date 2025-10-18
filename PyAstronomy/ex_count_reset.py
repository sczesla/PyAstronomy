from pathlib import Path

fns = list(Path(".").rglob("ex_*.ipynb"))

print("Files matching criteria:")
for i, fn in enumerate(fns,1):
    print(i, fn)
print()


for fn in fns:
    
    ec = 1
    ll = open(str(fn)).readlines()
    for i, l in enumerate(ll):
        if l.find('"execution_count":') != -1:
            ll[i] = l.split(":")[0] + ": %d,\n" % ec
            ec += 1

    open(str(fn), 'wt').writelines(ll)


