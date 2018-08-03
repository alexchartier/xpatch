import pickle
import glob

fnames = glob.glob('data/swarm/proc_lp/alex/*.pkl')
out_dir = 'data/pickle2/'
for fn in fnames:
    with open(fn, 'rb') as f:
        vals = pickle.load(f)
    out_fn = out_dir + fn.split('/')[-1]
    with open(out_fn, 'wb') as f:
        pickle.dump(vals, f, protocol=2)


