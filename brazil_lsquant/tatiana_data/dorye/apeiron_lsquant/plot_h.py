#!/usr/bin/env python3
import os, math, argparse
import numpy as np
import matplotlib.pyplot as plt

def read_wannier_hr_incomplete(fname):
    """
    Reads an incomplete Wannier90 hr.dat file:
      1) comment
      2) N = # orbitals
      3) M = # cells
      4) M degeneracy ints (skip)
      5) hopping lines: Rx Ry Rz i j Re Im
    Returns:
      N       : int
      Rs      : list of (rx,ry,rz)
      Present : bool array [i,j,R] where True if that entry was present
    """
    with open(fname, 'r') as f:
        _ = f.readline()                # comment
        N = int(f.readline().split()[0])
        M = int(f.readline().split()[0])
        for _ in range(math.ceil(M/15)):
            f.readline()                # skip degeneracies

        raw = []
        for line in f:
            parts = line.split()
            if len(parts) < 7: continue
            rx, ry, rz = map(int, parts[0:3])
            i,  j      = map(int, parts[3:5])
            raw.append((rx, ry, rz, i-1, j-1))

    if not raw:
        raise RuntimeError(f"No entries found in {fname!r}")

    # unique R‐vectors, sorted by ||R||
    Rs = sorted({(rx,ry,rz) for rx,ry,rz,_,_ in raw},
                key=lambda R: np.linalg.norm(R))
    R_to_idx = {R: k for k, R in enumerate(Rs)}

    # mark presence
    Present = np.zeros((N, N, len(Rs)), dtype=bool)
    for rx, ry, rz, i, j in raw:
        Present[i, j, R_to_idx[(rx, ry, rz)]] = True

    return N, Rs, Present

def main():
    p = argparse.ArgumentParser()
    p.add_argument("hr_file", help="Path to wannier90_hr.dat")
    args = p.parse_args()

    if not os.path.isfile(args.hr_file):
        print(f"Error: file not found: {args.hr_file!r}")
        exit(1)

    N, Rs, Present = read_wannier_hr_incomplete(args.hr_file)
    n_cells = len(Rs)

    # grid layout
    n_cols = int(math.ceil(math.sqrt(n_cells)))
    n_rows = int(math.ceil(n_cells / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3*n_cols, 3*n_rows),
                             constrained_layout=True,
                             squeeze=False)

    for idx, R in enumerate(Rs):
        r, c = divmod(idx, n_cols)
        ax   = axes[r][c]

        # show presence: True -> white (1), False -> black (0)
        img = Present[..., idx].astype(int)
        ax.imshow(img,
                  origin='lower',
                  cmap='gray',
                  interpolation='none',
                  aspect='equal',
                  vmin=0, vmax=1)

        ax.set_title(f"R = {R}", fontsize=9)
        ax.set_xlabel("j", fontsize=8)
        ax.set_ylabel("i", fontsize=8)

    # remove any empty subplots
    for idx in range(n_cells, n_rows*n_cols):
        r, c = divmod(idx, n_cols)
        fig.delaxes(axes[r][c])

    plt.suptitle("Presence map of H₍i,j₎(R): white=present, black=missing", fontsize=12)
    plt.savefig('test.png')

if __name__ == "__main__":
    main()
