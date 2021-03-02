import os


def maybe_gunzip(path):
    out = path[0:-3]
    if os.path.isfile(path) and path.endswith('.gz') and not os.path.isfile(out):
        os.system(f'gzcat {path} > {out}')
