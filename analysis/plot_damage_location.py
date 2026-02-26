"""
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  https://www.apache.org/licenses/LICENSE-2.0
"""

import numpy as np
import json
import sys
import matplotlib.pyplot as plt
import matplotlib.axes as axes
from typing import TypeAlias

# The atom name in a compound in a residue in a chain with damage type
DamageInfo: TypeAlias = tuple[str, int, str, str, int]
ELLIPSOID_COLORS = {'deoxyribose': '#BB5500', 'phosphate': "#C5C744", 'base': "#2BBC2B"}

def plot_damage_location(json_filename: str,
                         damage_filename: str):
    structure: json = read_json(json_filename)
    damages: list[DamageInfo] = read_damage(damage_filename)

    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')

    plot_structure(structure, ax)
    plot_damage(damages, structure, ax)

    ax.set_xlim(-30, 30)
    ax.set_ylim(-30, 30)
    ax.set_xlabel('x [Å]')
    ax.set_ylabel('z [Å]')
    ax.grid(True)
    plt.show()

def read_json(json_filename: str):
    with open(json_filename, mode="r") as f:
        structure = json.load(f)
    return structure

def read_damage(damage_filename: str):
    damages = []
    # skip headers
    with open(damage_filename, mode="r") as f:
        lines = f.readlines()[2:]
    # get each line
    for line in lines:
        line = line.rstrip('\n')
        if not line:
            continue
        parts = line.split()
        # chain_id, residue_id, compound_name, atom_name, damage_type
        damages.append((parts[2], int(parts[3]),
                        parts[4], parts[5], int(parts[6])))
    return damages

def plot_structure(structure: json, ax: axes.Axes):
    structure_center: np.ndarray = get_center(structure)
    
    for chain_info in structure["chains"]:
        for residue_info in chain_info["residues"]:
            for compound_info in residue_info["compounds"]:
                compound_name = compound_info["name"]
                # Get ellipsoidal parameters
                E = np.array(compound_info["ellipsoid"]["E"])
                center = np.array(compound_info["ellipsoid"]["center"]) - \
                    np.array([structure_center[0], structure_center[1], structure_center[2]])
                # projection of the ellipsoid to zx plane
                B, c_xz = compute_projection_quadratic(E, center)
                plot_ellipse_from_matrix(B, c_xz, ax, color=ELLIPSOID_COLORS[compound_name])

def plot_damage(damages: list[DamageInfo],
                structure: json, ax: axes.Axes):
    structure_center: np.ndarray = get_center(structure)
    hydrogen_name_dict = {'H1': 'H1\'',
                        'H2a': 'H2\'',
                        'H2b': 'H2\'\'',
                        'H3': 'H3\'',
                        'H4': 'H4\'',
                        'H5a': 'H5\'',
                        'H5b': 'H5\'\''}

    for damage in damages:
        chain_id = damage[0]
        residue_id = damage[1]
        compound_name = damage[2]
        atom_name = damage[3]
        damage_type = damage[4]
        
        if damage_type == 1: # direct damage
            atom_info = get_atom_info(chain_id, residue_id, compound_name, atom_name, structure)
            if atom_info is not None:
                atom_position = np.array(atom_info["position"]) - \
                    np.array([structure_center[0], structure_center[1], structure_center[2]])
                ax.plot(atom_position[0], atom_position[2], 'ro', markersize=1)

        elif damage_type == 2: # indirect damage
            hydrogen_name = hydrogen_name_dict[atom_name]
            atom_info = get_atom_info(chain_id, residue_id, compound_name, hydrogen_name, structure)
            if atom_info is not None:
                atom_position = np.array(atom_info["position"]) - \
                    np.array([structure_center[0], structure_center[1], structure_center[2]])
                ax.plot(atom_position[0], atom_position[2], 'bo', markersize=1)

def get_atom_info(chain_id: str, residue_id: int, compound_name: str,
                  atom_name: str, structure: json):
    for chain_info in structure["chains"]:
        if chain_info["id"] != chain_id:
            continue
        for residue_info in chain_info["residues"]:
            if residue_info["id"] != residue_id:
                continue
            for compound_info in residue_info["compounds"]:
                if compound_info["name"] != compound_name:
                    continue
                for atom_info in compound_info["atoms"]:
                    if atom_info["name"] == atom_name:
                        return atom_info
    return None

def get_center(structure: json):
    x_vals = []
    y_vals = []
    z_vals = []

    for chain_info in structure["chains"]:
        for residue_info in chain_info["residues"]:
            for compound_info in residue_info["compounds"]:
                center = np.array(compound_info["ellipsoid"]["center"])
                x_vals.append(center[0])
                y_vals.append(center[1])
                z_vals.append(center[2])
    x_max, x_min = max(x_vals), min(x_vals)
    y_max, y_min = max(y_vals), min(y_vals)
    z_max, z_min = max(z_vals), min(z_vals)

    # center of gravity
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    z_center = (z_max + z_min) / 2

    return np.array([x_center, y_center, z_center])

def compute_projection_quadratic(E: np.ndarray, c: float) -> np.ndarray:
    """
    Calculates the projected 2-dim matrix B and center c_xz into y axis direction
    from the 3-dim matrix E representing the ellipsoidal shepe and rotation.
    """
    # Block separation
    A = E[np.ix_([0, 2], [0, 2])]   # row 0,2 col 0,2
    b = E[np.ix_([0, 2], [1])]      # row 0,2 col 1
    e22 = E[1, 1]

    # projected matrix B, center
    B = A - (b @ b.T) / e22
    c_xz = np.array([c[0], c[2]])

    return B, c_xz

def plot_ellipse_from_matrix(
    B: np.ndarray,
    c: np.ndarray,
    ax: axes.Axes,
    num_points: int = 200,
    level: float = 1.0,
    color: str = 'blue',
):

    eigvals, eigvecs = np.linalg.eigh(B)
    axes_lengths = np.sqrt(level / eigvals)
    rotation = eigvecs

    t = np.linspace(0, 2 * np.pi, num_points)
    circle = np.vstack([np.cos(t), np.sin(t)])
    ellipse_points = rotation @ (np.diag(axes_lengths) @ circle)
    ellipse_points[0, :] += c[0]
    ellipse_points[1, :] += c[1]

    ax.fill(ellipse_points[0, :], ellipse_points[1, :], color=color, alpha=0.1)
    ax.plot(ellipse_points[0, :], ellipse_points[1, :], color=color, linewidth=1.0)

    return ax

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Please specify the input structure and damage file (.txt).")
        print("Usage example: python3 print_damage_summary.py \n \
                                      data/json/nuc-SASA.json\n \
                                      data/strand_breaks/StrandBreakInfo.txt")

    json_filename = sys.argv[1]
    damage_filename = sys.argv[2]

    plot_damage_location(json_filename, damage_filename)