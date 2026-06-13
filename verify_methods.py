import re
import subprocess

with open('vignettes/manuscript/chapter1.qmd') as f:
    new = f.read()

result = subprocess.run(
    ['git', 'show', 'HEAD:vignettes/manuscript/chapter1.qmd'],
    capture_output=True, text=True
)
old = result.stdout

def normalize(text):
    return re.sub(r'\s+', ' ', text).strip()

if normalize(old) == normalize(new):
    print('IDENTICAL after whitespace normalization')
else:
    o = normalize(old).split()
    n = normalize(new).split()
    for i, (a, b) in enumerate(zip(o, n)):
        if a != b:
            print(f'First diff at word {i}: old={repr(a)}, new={repr(b)}')
            print('Context old:', ' '.join(o[max(0,i-5):i+5]))
            print('Context new:', ' '.join(n[max(0,i-5):i+5]))
            break
    else:
        print(f'Length differs: old={len(o)}, new={len(n)}')
