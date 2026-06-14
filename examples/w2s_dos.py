#!/usr/bin/env python3
"""Load a wannier2sparse CSR file and plot its spectral density (DOS) via KPM,
or the sorted eigenvalue spectrum. Only needs numpy/scipy/matplotlib.

Usage:
    python3 w2s_dos.py PATH.HAM.CSR [--title T] [--out fig.png]
                       [--mode dos|spectrum] [--moments 2048] [--vectors 20]

The CSR text format written by wannier2sparse is:
    line 1: dim nnz
    line 2: 2*nnz reals  -> (real, imag) of each stored value, row-major
    line 3: nnz ints     -> column index of each value
    line 4: dim+1 ints   -> row pointers
"""
import sys, argparse, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def load_csr(path):
    with open(path) as f:
        dim, nnz = map(int, f.readline().split())
        v   = np.fromstring(f.readline(), sep=' '); val = v[0::2] + 1j*v[1::2]
        col = np.fromstring(f.readline(), sep=' ', dtype=np.int64)
        rp  = np.fromstring(f.readline(), sep=' ', dtype=np.int64)
    H = sp.csr_matrix((val, col, rp), shape=(dim, dim))
    return H

def kpm_dos(H, M=2048, R=20, npts=2000, pad=1.05, seed=0):
    """Kernel Polynomial Method DOS: Chebyshev moments by stochastic trace,
    Jackson kernel. Uses only sparse mat-vec products (scales to large dim)."""
    n = H.shape[0]
    emax = float(eigsh(H, k=1, which='LA', return_eigenvectors=False)[0])
    emin = float(eigsh(H, k=1, which='SA', return_eigenvectors=False)[0])
    a = (emax - emin) / 2 * pad
    b = (emax + emin) / 2
    Hs = (H - b * sp.identity(n, format='csr')) / a       # rescaled to [-1,1]
    mu = np.zeros(M)
    rng = np.random.default_rng(seed)
    for _ in range(R):
        v = np.exp(2j*np.pi*rng.random(n)); v /= np.linalg.norm(v)
        t0 = v; t1 = Hs @ t0
        mu[0] += np.vdot(v, t0).real
        mu[1] += np.vdot(v, t1).real
        for m in range(2, M):
            t2 = 2*(Hs @ t1) - t0
            mu[m] += np.vdot(v, t2).real
            t0, t1 = t1, t2
    mu /= R
    m = np.arange(M)                                       # Jackson kernel
    g = ((M-m+1)*np.cos(np.pi*m/(M+1)) + np.sin(np.pi*m/(M+1))/np.tan(np.pi/(M+1)))/(M+1)
    mu *= g
    x = np.linspace(-1, 1, npts + 2)[1:-1]  # open interval: avoid 1/sqrt(1-x^2) blow-up at +-1
    T = np.cos(np.outer(np.arccos(np.clip(x, -1, 1)), m))
    rho = (mu[0] + 2*np.sum(mu[1:]*T[:, 1:], axis=1)) / (np.pi*np.sqrt(1 - x**2))
    return a*x + b, rho/(a*n)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('csr'); p.add_argument('--title', default='')
    p.add_argument('--out', default='dos.png'); p.add_argument('--mode', default='dos',
                   choices=['dos', 'spectrum'])
    p.add_argument('--moments', type=int, default=2048)
    p.add_argument('--vectors', type=int, default=20)
    args = p.parse_args()

    H = load_csr(args.csr)
    herr = abs(H - H.getH()).max() if H.nnz else 0.0
    print("loaded %s  dim=%d  nnz=%d  hermiticity_err=%.2e" % (args.csr, H.shape[0], H.nnz, herr))

    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    if args.mode == 'spectrum':                            # exact, for SMALL supercells
        E = np.linalg.eigvalsh(H.toarray())
        plt.plot(np.linspace(0, 1, E.size), E, '.', ms=2)
        plt.xlabel('state index (normalized)'); plt.ylabel('E')
    else:
        E, rho = kpm_dos(H, M=args.moments, R=args.vectors)
        plt.fill_between(E, rho, alpha=.4); plt.plot(E, rho)
        plt.xlabel('E'); plt.ylabel('DOS')
    plt.title(args.title or args.csr); plt.tight_layout(); plt.savefig(args.out, dpi=120)
    print("saved", args.out)

if __name__ == '__main__':
    main()
