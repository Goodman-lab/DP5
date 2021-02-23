import gzip
import pickle


atomic_reps = pickle.load(open("frag_reps.p","rb"))

with gzip.open("frag_reps.gz", "wb") as f:
    pickle.dump(atomic_reps, f)


