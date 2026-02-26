"""
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
"""

import sys
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

candidates = [
    "Hiragino Sans",         # mac
    "Hiragino Kaku Gothic Pro",
    "IPAexGothic",           # Linux
    "IPAGothic", 
    "Noto Sans CJK JP",      # CJK
    "Arial Unicode MS"       # Unicode
]

available = {f.name for f in fm.fontManager.ttflist}

for name in candidates:
    if name in available:
        plt.rcParams["font.family"] = name
        print("Font Used:", name)
        break

# G value with RMS dictionary at each time of each species
SpeciesInfo = dict[str, list[tuple[float, float, float]]]

def main():
    args = sys.argv
    if len(args) == 1:
        print("Please provide a valid .txt filename.")
        print("Usage Example: python3 plot_G_value.py data/speceis/SpeciesInfo.txt")
        exit(1)
    
    plot_G_value(args[1])

def plot_G_value(filename: str):
    species_info = get_species_info(filename)
    plot_species_G_value(species_info)
    plot_damage_G_value(species_info)

def get_species_info(filename: str):
    species_info: SpeciesInfo = {}
    with open(filename, mode="r") as f:
        for line in f:
            # if a line starts from "#", skip it
            if line[0] == "#":
                continue
            attrs = line.split()
            time_ps = float(attrs[0])
            g_value = float(attrs[1])
            rms = float(attrs[2])
            mole_name = attrs[4]

            if mole_name not in species_info:
                species_info[mole_name] = []
            species_info[mole_name].append((time_ps, g_value, rms))

    return species_info

def plot_species_G_value(species_info: SpeciesInfo):
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")

    plot_species = [
        "e_aq^-1",
        "°OH^0",
        "H^0",
        "H2O2^0",
        "HO_2°^0",
        "O_2^-1"
    ]
    markers = [
        ('o', False),
        ('o', True),
        ('^', False),
        ('^', True),
        ('s', False),
        ('s', True)
    ]

    for idx, species in enumerate(plot_species):
        info = species_info.get(species)
        shape, fill = markers[idx]
        if fill:
            back_col = 'black'
        else:
            back_col = 'none'
        if info is not None:
            time_ps_list = [time_ps for time_ps,_,_ in info]
            g_value_list = [g_value for _,g_value,_ in info]
            rms_list = [rms for _,_,rms in info]

            plt.errorbar(
                        time_ps_list,
                        g_value_list,
                        rms_list,
                        marker=shape,
                        linestyle='-',
                        color='black',
                        ecolor='black',
                        markeredgecolor='black',
                        markerfacecolor=back_col,
                        markersize=3,
                        linewidth=0.5,
                        capsize=1
                        )

    ax.set_xlim(1e0, 1e7)
    ax.set_xlabel("Time [ps]", fontsize=16, loc='right',
                  labelpad=15)
    ax.set_ylabel("G Value [molecules/100 eV]", fontsize=16, loc='top',
                  labelpad=15)
    ax.tick_params(labelsize=18)
    ax.minorticks_on()
    plt.show()
    plt.close(fig)

def plot_damage_G_value(species_info: SpeciesInfo):

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")

    plot_species = [
        "e_aq^-1",
        "°OH^0",
        "H^0",
        "H2O2^0",
        "HO_2°^0",
        "O_2^-1"
    ]
    markers = [
        ('o', False),
        ('o', True),
        ('^', False),
        ('^', True),
        ('s', False),
        ('s', True),
        ('*', False)
    ]

    plot_species = [
        "DNA_DamagedDeoxyribose_H1_OH^0",
        "DNA_DamagedDeoxyribose_H2a_OH^0",
        "DNA_DamagedDeoxyribose_H2b_OH^0",
        "DNA_DamagedDeoxyribose_H3_OH^0",
        "DNA_DamagedDeoxyribose_H4_OH^0",
        "DNA_DamagedDeoxyribose_H5a_OH^0",
        "DNA_DamagedDeoxyribose_H5b_OH^0"
    ]

    for idx, species in enumerate(plot_species):
        info = species_info.get(species)
        shape, fill = markers[idx]
        if fill:
            back_col = 'black'
        else:
            back_col = 'none'
        if info is not None:
            time_ps_list = [time_ps for time_ps,_,_ in info]
            g_value_list = [g_value for _,g_value,_ in info]
            rms_list = [rms for _,_,rms in info]

            plt.errorbar(
                        time_ps_list,
                        g_value_list,
                        rms_list,
                        marker=shape,
                        linestyle='-',
                        color='black',
                        ecolor='black',
                        markeredgecolor='black',
                        markerfacecolor=back_col,
                        markersize=3,
                        linewidth=0.5,
                        capsize=1
                        )

    ax.set_xlim(1e0, 1e7)
    ax.set_xlabel("Time [ps]", fontsize=16, loc='right',
                  labelpad=15)
    ax.set_ylabel("G Value [molecules/100 eV]", fontsize=16, loc='top',
                  labelpad=5)
    ax.tick_params(labelsize=18)
    ax.minorticks_on()
    fig.subplots_adjust(
        left=0.15,
        right=0.95,
        bottom=0.15,
        top=0.95
        )
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    main()