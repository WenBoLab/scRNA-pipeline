import os
for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"]:
    os.environ[key] = "1"
os.environ["MPLBACKEND"] = "Agg"

