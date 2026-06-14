#!/usr/bin/env python3
"""Generate canonical tight-binding models as Wannier90 inputs for wannier2sparse.
Writes <model>/<model>{_hr.dat,.uc,.xyz} for: chain1d, graphene, cubic, haldane.
Each model has a known analytic spectrum (see examples/README.md)."""
import os, math

def write(model, num_wann, rpts, hoppings, uc, xyz):
    os.makedirs(model, exist_ok=True)
    with open(f"{model}/{model}_hr.dat", "w") as f:
        f.write(f"{model} (generated)\n{num_wann}\n{len(rpts)}\n")
        f.write(" ".join("1" for _ in rpts) + "\n")          # ndegen = 1 everywhere
        for (R, i, j, re, im) in hoppings:
            f.write("%d %d %d  %d %d  %.10f %.10f\n" % (R[0], R[1], R[2], i, j, re, im))
    open(f"{model}/{model}.uc", "w").write(uc)
    open(f"{model}/{model}.xyz", "w").write(xyz)
    print("wrote", model)

# --- 1D chain: 1 orbital, NN t=-1.  E(k)=2t cos k -> band [-2,2] -----------
def chain1d():
    h = [((-1,0,0),1,1,-1,0), ((0,0,0),1,1,0,0), ((1,0,0),1,1,-1,0)]
    write("chain1d", 1, {(-1,0,0),(0,0,0),(1,0,0)}, h,
          "1 0 0\n0 1 0\n0 0 1\n", "1\nA 0 0 0\n")

# --- graphene: 2 sublattices, NN t=-1.  Dirac point at E=0, edges +-3 -------
def graphene():
    t=-1
    h = [((0,0,0),1,2,t,0), ((0,0,0),2,1,t,0),
         ((-1,0,0),1,2,t,0), ((1,0,0),2,1,t,0),
         ((0,-1,0),1,2,t,0), ((0,1,0),2,1,t,0)]
    write("graphene", 2, {(0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)}, h,
          "1.5 0.8660254 0\n1.5 -0.8660254 0\n0 0 10\n", "2\nCA 0 0 0\nCB 1 0 0\n")

# --- simple cubic: 1 orbital, NN t=-1 in x,y,z.  edges +-6, van Hove --------
def cubic():
    t=-1; R=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    h=[((0,0,0),1,1,0,0)] + [(r,1,1,t,0) for r in R]
    write("cubic", 1, set([(0,0,0)]+R), h,
          "1 0 0\n0 1 0\n0 0 1\n", "1\nA 0 0 0\n")

# --- Haldane: graphene + complex NNN t2 e^{i phi}.  Opens a gap at E=0 ------
def haldane(t1=-1.0, t2=0.15, phi=math.pi/2):
    c, s = t2*math.cos(phi), t2*math.sin(phi)
    h = [((0,0,0),1,2,t1,0), ((0,0,0),2,1,t1,0),
         ((-1,0,0),1,2,t1,0), ((1,0,0),2,1,t1,0),
         ((0,-1,0),1,2,t1,0), ((0,1,0),2,1,t1,0)]
    for R in [(1,0,0),(0,-1,0),(-1,1,0)]:                    # A: +phi, B: -phi
        h += [(R,1,1,c,s), (R,2,2,c,-s)]
    for R in [(-1,0,0),(0,1,0),(1,-1,0)]:                    # conjugates
        h += [(R,1,1,c,-s), (R,2,2,c,s)]
    rpts = set(t[0] for t in h)
    write("haldane", 2, rpts, h,
          "1.5 0.8660254 0\n1.5 -0.8660254 0\n0 0 10\n", "2\nCA 0 0 0\nCB 1 0 0\n")

# --- magnetized 1D chain: 1 site, 2 spin orbitals, NN t=-1, on-site exchange -J_ex sigma_z.
#     sup/sdw convention: orbital labels carry _s+_ (up) / _s-_ (down) so the spin
#     operators (S_x,S_y,S_z) are built from the spin doubling. Orbital 1 = up, 2 = down;
#     bands E_sigma(k) = -2 cos(k) -/+ J_ex  ->  spectrum [-2-J_ex, 2+J_ex]. -----------
def magnetic_chain(Jex=0.1):
    h = [((0,0,0),1,1,-Jex,0), ((0,0,0),2,2,Jex,0),
         ((-1,0,0),1,1,-1,0),  ((-1,0,0),2,2,-1,0),
         ((1,0,0),1,1,-1,0),   ((1,0,0),2,2,-1,0)]
    write("chain1d_mag", 2, {(-1,0,0),(0,0,0),(1,0,0)}, h,
          "1 0 0\n0 1 0\n0 0 1\n", "2\nA_s+_ 0 0 0\nA_s-_ 0 0 0\n")

if __name__ == "__main__":
    chain1d(); graphene(); cubic(); haldane(); magnetic_chain()
