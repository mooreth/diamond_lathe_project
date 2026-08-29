#!/usr/bin/env python3
"""
convert_diamond_fixture.py
---------------------------------------------------------------
Converts the SAMSON-exported PDB of the diamond lathe fixture
(rail, sliding tool, two rotating grippers) into a LAMMPS data
file, keeping each PDB chain as a separate LAMMPS molecule-ID so
that the pieces can be grouped and driven independently in the
LAMMPS input script (fix move / fix rigid).

Chain -> part mapping (inferred geometrically, see README block
printed at the end -- please verify visually in OVITO by
coloring on "Molecule Identifier" before trusting it blindly):

    A, E        tool holder + cutting-tool tip  (slides along rail, Y)
    B, C, D     rail beam + its two end mounts   (stationary, or one
                bulk infeed move in X/Z, then held fixed)
    F, G, H     right gripper fingers
    I           right gripper rotating hub/disc
    J           right collet nose (clamps workpiece)
    K, L, M     left gripper fingers
    N           left gripper rotating hub/disc
    O           left collet nose (clamps workpiece)

Only C and H atoms are present in the file (hydrogen-terminated
diamond parts). The soft-metal workpiece itself is NOT in this
file -- it is created separately inside the LAMMPS script as an
EAM lattice (Cu, using CuAlW.txt) threaded along the rotation
axis, in the gap between the two collets.

Usage:
    python3 convert_diamond_fixture.py lathe_diamond_parts.pdb diamond_fixture.data
"""

import sys
import numpy as np

MASS = {"C": 12.011, "H": 1.008}
TYPE = {"C": 1, "H": 2}


def parse_pdb(path):
    chains = {}          # chain_letter -> list[(x,y,z,elem)]
    order = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip() or line[12:14].strip()[0]
                if chain not in chains:
                    chains[chain] = []
                    order.append(chain)
                chains[chain].append((x, y, z, elem))
    return chains, sorted(order)


def main():
    if len(sys.argv) != 3:
        print("usage: convert_diamond_fixture.py input.pdb output.data")
        sys.exit(1)

    in_pdb, out_data = sys.argv[1], sys.argv[2]
    chains, order = parse_pdb(in_pdb)

    mol_id = {c: i + 1 for i, c in enumerate(order)}  # A=1, B=2, ... O=15

    all_pts = []
    atoms = []   # (chain, x, y, z, elem)
    for c in order:
        for (x, y, z, e) in chains[c]:
            atoms.append((c, x, y, z, e))
            all_pts.append((x, y, z))
    all_pts = np.array(all_pts)

    # Padding around the fixture. Rotation of the grippers stays within
    # their own original radial envelope, so we mainly need padding for
    # (a) the tool's small radial infeed and (b) a safety margin against
    # atoms sitting exactly on a fixed boundary.
    pad = 25.0
    lo = all_pts.min(axis=0) - pad
    hi = all_pts.max(axis=0) + pad

    with open(out_data, "w") as f:
        f.write("LAMMPS data file: diamond lathe fixture "
                "(rail + tool + two rotating grippers)\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write("3 atom types\n\n")   # 1=C 2=H reserved, 3=Cu added later in LAMMPS
        f.write(f"{lo[0]:.4f} {hi[0]:.4f} xlo xhi\n")
        f.write(f"{lo[1]:.4f} {hi[1]:.4f} ylo yhi\n")
        f.write(f"{lo[2]:.4f} {hi[2]:.4f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {MASS['C']:.4f}  # C (diamond)\n")
        f.write(f"2 {MASS['H']:.4f}  # H (diamond surface termination)\n")
        f.write(f"3 63.5460          # Cu (workpiece, added later)\n\n")

        f.write("Atoms # full\n\n")
        for i, (c, x, y, z, e) in enumerate(atoms, start=1):
            t = TYPE[e]
            m = mol_id[c]
            f.write(f"{i} {m} {t} 0.0 {x:.4f} {y:.4f} {z:.4f}\n")

    # ---- report geometry so the LAMMPS script's numbers can be checked ----
    def bbox(cs):
        pts = np.array([(x, y, z) for c, x, y, z, e in atoms if c in cs])
        return pts.min(axis=0), pts.max(axis=0)

    print(f"Wrote {out_data}: {len(atoms)} atoms, "
          f"{len(order)} chains -> molecules 1..{len(order)}")
    print(f"Chain -> molecule-ID map: "
          f"{', '.join(f'{c}={mol_id[c]}' for c in order)}")
    print(f"Box (with {pad} A padding): "
          f"x[{lo[0]:.1f},{hi[0]:.1f}] "
          f"y[{lo[1]:.1f},{hi[1]:.1f}] "
          f"z[{lo[2]:.1f},{hi[2]:.1f}]")

    for label, cs in [("rail (B,C,D)", "BCD"),
                       ("tool (A,E)", "AE"),
                       ("gripper_right (F,G,H,I,J)", "FGHIJ"),
                       ("gripper_left (K,L,M,N,O)", "KLMNO")]:
        mn, mx = bbox(set(cs))
        print(f"  {label:28s} bbox lo={mn.round(1)} hi={mx.round(1)}")


if __name__ == "__main__":
    main()
