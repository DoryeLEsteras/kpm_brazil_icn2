import os
import numpy as np
import matplotlib.pyplot as plt
import pybinding as pb
from pybinding import Series
from math import sqrt

# ─── user‑adjustable constants ─────────────────────────
ORBITALS_PER_ATOM = 18
RYDBERG_TO_EV     = 1.0       # set to 13.605693 if _hr.dat is in Rydberg
STEP_K            = 0.025     # k‑path resolution
DOS_BROADENING    = 0.01      # in eV
DOS_NPTS          = 10_000    # points in the energy grid
KPOINTS_GRID      = 400        # ← define your Monkhorst–Pack grid size here
ENERGY_WINDOW     = (-20, +20) # (emin, emax) around E_F, in eV
DOS_THRESHOLD     = 1e-3      # to detect band edges in DOS
# ───────────────────────────────────────────────────────────────────

def parse_wannier_data(path, prefix):
    """Read <prefix>_lat.dat and <prefix>_hr.dat. 
    Returns: lat_vec(3×3), positions list, Ef (float), ham_r dict."""
    # 1) read + clean lines
    latfile = os.path.join(path, f"{prefix}_lat.dat")
    with open(latfile) as f:
        lines = [L.strip() for L in f if L.strip()]
    # 2) lattice vectors
    i0, i1 = lines.index("begin unit_cell")+1, lines.index("end unit_cell")
    ucell = lines[i0:i1]
    if len(ucell)!=3:
        raise RuntimeError("need 3 lattice vectors")
    lat_vec = np.zeros((3,3))
    for i, L in enumerate(ucell):
        parts = [float(x) for x in L.split(',') if x.strip()]
        lat_vec[i] = parts
    # 3) atomic positions
    p0, p1 = lines.index("begin position")+1, lines.index("end position")
    pos_lines = lines[p0:p1]
    positions = []
    for L in pos_lines:
        sl = [x.strip() for x in L.split(',') if x.strip()]
        label, x,y,z = sl[0], *map(float, sl[1:])
        positions.append((label, x,y,z))
    # 4) Fermi energy
    f0, f1 = lines.index("begin fermi")+1, lines.index("end fermi")
    Ef_line = lines[f0:f1]
    if len(Ef_line)!=1:
        raise RuntimeError("need exactly one Fermi energy")
    Ef = float(Ef_line[0])
    # 5) hoppings: hr.dat
    hrfile = os.path.join(path, f"{prefix}_hr.dat")
    num_orb = len(positions)*ORBITALS_PER_ATOM
    ham_r = {}
    with open(hrfile) as f:
        for L in f:
            if not L.strip() or L.startswith('#'):
                continue
            sp = L.split()
            R = tuple(map(int, sp[0:3]))
            i, j = int(sp[3])-1, int(sp[4])-1
            val = (float(sp[5]) + 1j*float(sp[6])) * RYDBERG_TO_EV
            if R not in ham_r:
                ham_r[R] = np.zeros((num_orb,num_orb), dtype=complex)
            ham_r[R][i,j] = val
    return lat_vec, positions, Ef, ham_r

def build_model_2d(lat_vec, positions, ham_r):
    """Constructs the 2D pb.Model and returns (model2d, b1, b2)."""
    num_orb = len(positions)*ORBITALS_PER_ATOM
    # onsite energies from R=(0,0,0)
    onsite = np.real(np.diag(ham_r[(0,0,0)]))
    # 3D lattice
    lat3d = pb.Lattice(a1=lat_vec[0], a2=lat_vec[1], a3=lat_vec[2])
    # map each orbital → a sublattice name
    orbital_labels = {}
    for a_idx, (lab, x,y,z) in enumerate(positions):
        for orb in range(ORBITALS_PER_ATOM):
            idx = a_idx*ORBITALS_PER_ATOM + orb
            name = f"{lab}{orb:02d}"
            orbital_labels[idx] = name
            lat3d.add_one_sublattice(name, [x,y,z], onsite[idx])
    # hoppings (only add each ±R pair once)
    added = []
    for R, H in ham_r.items():
        if R in added or tuple(-np.array(R)) in added:
            continue
        for i in range(num_orb):
            for j in range(num_orb):
                if R==(0,0,0) and i>=j:
                    continue
                val = H[j,i]
                if val != 0:
                    lat3d.add_one_hopping(R, orbital_labels[i], orbital_labels[j], val)
        added.append(R)
    # project to 2D
    a12d = lat3d.vectors[0][:2]
    a22d = lat3d.vectors[1][:2]
    lat2d = pb.Lattice(a1=a12d, a2=a22d)
    # copy sublattices
    for name, sub in lat3d.sublattices.items():
        lat2d.add_one_sublattice(name, sub.position[:2], np.real(sub.energy[0][0]))
    # copy hoppings
    for _, hop in lat3d.hoppings.items():
        for term in hop.terms:
            cell = term.relative_index[:2]
            lat2d.add_one_hopping(
                cell,
                orbital_labels[term.from_id],
                orbital_labels[term.to_id],
                hop.energy[0][0]
            )
    model2d = pb.Model(lat2d, pb.translational_symmetry(a1=True, a2=True))
    # reciprocal vectors (full 3D)
    a1,a2,a3 = lat_vec
    b1 = 2*np.pi * np.cross(a2,a3) / np.dot(a1, np.cross(a2,a3))
    b2 = 2*np.pi * np.cross(a3,a1) / np.dot(a2, np.cross(a3,a1))
    return model2d, b1, b2

def calculate_bands(model2d, b2, Ef, step=STEP_K):
    """Returns (bands, solver) for the 2D Γ–M–K–Γ path."""
    G  = [0, 0]
    M  = [0, 0.5*b2[1]]
    K  = [0.5*b2[1]/sqrt(3), 0.5*b2[1]]
    solver = pb.solver.lapack(model2d)
    bands  = solver.calc_bands(G, M, K, G, step=step)
    print(bands)
    return bands, solver

def plot_bands(bands, Ef, filename, window=ENERGY_WINDOW):
    x = bands.k_path.as_1d()
    y = bands.energy - Ef
    plt.figure()
    for band in y.T:
        plt.plot(x, band, lw=2)
    # high‑symmetry lines
    xlocs = x[bands.k_path.point_indices]
    plt.xticks(xlocs, [r"$\Gamma$", "M", "K", r"$\Gamma$"])
    for xv in xlocs:
        plt.axvline(xv, color='gray', linestyle='--', linewidth=1)
    plt.ylim(window)
    plt.xlim(x[0], x[-1])
    yt = np.arange(window[0], window[1]+1, 1.0)
    plt.yticks(yt)
    plt.ylabel("E − E$_F$ (eV)")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)

def calculate_dos(solver, b1, b2, Ef,
                  window=ENERGY_WINDOW,
                  npts=DOS_NPTS,
                  broadening=DOS_BROADENING,
                  nk=KPOINTS_GRID):
    """Manual DOS via k‑point summation on a Monkhorst–Pack nk×nk grid."""
    emin, emax = window
    # absolute energy grid
    energies = np.linspace(Ef + emin, Ef + emax, npts)
    dos_accum = np.zeros_like(energies)

    # fractional k‑grid
    k_fracs = np.linspace(0, 1, nk, endpoint=False)
    for kx in k_fracs:
        for ky in k_fracs:
            kvec = kx * b1[:2] + ky * b2[:2]
            solver.set_wave_vector([kvec[0], kvec[1]])
            eigs = solver.eigenvalues
            # Gaussian broadening
            dos_accum += np.exp(
                - (energies[:,None] - eigs[None,:])**2 / (2 * broadening**2)
            ).sum(axis=1)

    # normalize by BZ area and Gaussian prefactor
    area_bz      = abs(np.cross(b1[:2], b2[:2]))
    gauss_factor = broadening * np.sqrt(2*np.pi)
    dos_accum *= area_bz / (nk*nk * gauss_factor)

    # return a simple Series-like object
    return Series(variable=energies, data=dos_accum)

def plot_dos(dos, Ef, filename, window=ENERGY_WINDOW):
    E_rel = dos.variable - Ef
    plt.figure()
    plt.plot(E_rel, dos.data, lw=2)
    plt.axvline(0, color='gray', linestyle='--', linewidth=1)
    plt.xlim(window)
    xt = np.arange(window[0], window[1]+1, 1.0)
    plt.xticks(xt)
    plt.xlabel("E − E$_F$ (eV)")
    plt.ylabel("DOS (states/eV)")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)

def compute_band_gap(bands, Ef):
    E_rel = (bands.energy - Ef).ravel()
    vb_max = E_rel[E_rel < 0].max()
    cb_min = E_rel[E_rel > 0].min()
    return (cb_min - vb_max), vb_max, cb_min

def compute_dos_gap(dos, Ef, threshold=DOS_THRESHOLD):
    E_rel = dos.variable - Ef
    mask  = dos.data > threshold
    vb_e = E_rel[(mask) & (E_rel < 0)].max()
    cb_e = E_rel[(mask) & (E_rel > 0)].min()
    return (cb_e - vb_e), vb_e, cb_e


def main():
    # — user parameters —
    folder   = "./"
    seedname = "mos2"
    path, pre = folder, seedname

    # parse & build
    lat_vec, positions, Ef, ham_r = parse_wannier_data(path, pre)
    model2d, b1, b2              = build_model_2d(lat_vec, positions, ham_r)

    # bands
    bands, solver = calculate_bands(model2d, b2, Ef)
    plot_bands(bands, Ef, f"{seedname}_bands.png")
    bg, vb, cb   = compute_band_gap(bands, Ef)
    print(f"Band‐structure gap = {bg:.3f} eV (VB @ {vb:.3f}, CB @ {cb:.3f})")

    # DOS
    dos = calculate_dos(solver, b1, b2, Ef, nk=KPOINTS_GRID)
    plot_dos(dos, Ef, f"{seedname}_dos.png")
    dg, vbe, cbe = compute_dos_gap(dos, Ef)
    print(f"DOS gap            = {dg:.3f} eV (VB @ {vbe:.3f}, CB @ {cbe:.3f})")

if __name__ == "__main__":
    main()
    