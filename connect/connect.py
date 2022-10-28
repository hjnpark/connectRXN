#! /usr/bin/env python
import numpy as np
import sys, os, copy, itertools

from pyvis.network import Network
from calc_rmsd.permutation_invariant_rmsd import min_rmsd
from .molecule import Molecule, TopEqual, Elements
from collections import OrderedDict
import tarfile
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import colors

from xyz2mol import xyz2mol
from rdkit import Chem
from rdkit.Chem import Draw

au2kcal = 627.5094740630558
au2kj = 2625.4996394798254

def equal(m1, m2, threshold=0.05):
    """
    To see whether the two Molecule objects have the same geometry

    Parameters
    ----------
    m1, m2: Molecule Object

    threshold: float

    Return
    ----------
    True or False
    """
    top_equal = TopEqual(m1, m2)
    if not top_equal:
        return False
    minimum_rmsd = min_rmsd(m1, m2)[0]

    # To check whether this function works well, writes xyz file so their geometries can be manually investigated.
    #import random

    #rand_int = random.randint(0, 100)
    #if minimum_rmsd < threshold and rand_int == 50:
    #    M = m1 + m2
    #    M.align()
    #    if not os.path.exists("low_rmsd"):
    #        os.mkdir("low_rmsd")
    #    M.comms = [
    #        "Minimum RMSD between the two molecules is %f Angstrom" % minimum_rmsd,
    #        "",
    #    ]
    #    mol_num = (len(os.listdir("low_rmsd"))) + 1
    #    M.write("low_rmsd/%i.xyz" % mol_num)

    return minimum_rmsd < threshold


class BuildGraph(object):
    """
    Use networkx to build a graph
    """

    def __init__(self, rxns, charge):
        self.G = nx.Graph()
        self.rxns = rxns
        self.charge = charge
        self.unique_rxns = []
        self.lowest_E = 0.0
        self.highest_E = 0.0

    def unify(self):
        """Pick out same molecules with different molecule objects and unify them"""
        print("Detecting different molecule objects with same geometries and unifying them.")
        rxn_list = copy.deepcopy(list(self.rxns.values()))
        unique_pairs = list(itertools.combinations(enumerate(rxn_list), 2))
        skipping = []
        for (i, u_rxn1), (j, u_rxn2) in unique_pairs:
            M_list = [u_rxn1[0], u_rxn1[-1], u_rxn2[0], u_rxn2[-1]]
            M_index = [0, -1, 0, -1]
            TS1, TS2 = u_rxn1[1], u_rxn2[1]
            for (k, M1), (l, M2) in list(itertools.combinations(enumerate(M_list), 2)):
                if M1 not in skipping:
                    if M1.qm_energies[0] < self.lowest_E or M1.qm_energies[0] < self.lowest_E:
                        self.lowest_E = min(M1.qm_energies[0], M2.qm_energies[0])
                    if M1.qm_energies[0] < self.highest_E or M1.qm_energies[0] < self.highest_E:
                        self.highest_E = min(M1.qm_energies[0], M2.qm_energies[0])

                    if equal(M1, M2):
                        rxn_list[j][M_index[l]] = M1
                        skipping.append(M2)

            if TS1 not in skipping and equal(TS1, TS2):
                rxn_list[j][1] = TS1
                skipping.append(TS2)
        self.unique_rxns = rxn_list

    def build(self):
        temp_G = nx.Graph()
        self.unify()
        print("Building the graph...")
        for v in self.unique_rxns:
            for i, M in enumerate(v):
                if i % 2 != 1:
                    atoms = []
                    for atom in v[i].elem:
                        atoms.append(Elements.index(atom))
                    a = xyz2mol(
                        atoms, v[i].xyzs[0], charge=self.charge
                    )  # , allow_charged_fragments=False)
                    if len(a) != 0:
                        smiles = Chem.MolToSmiles(a[0])

                        if not os.path.exists("2D_images"):
                            os.mkdir("2D_images")
                        Draw.MolToFile(
                            a[0],
                            os.path.join(
                                "2D_images", "%s.png" % smiles.replace("/", "_")
                            ),
                        )

                    else:
                        print(
                            "WARNING: A node is missing SMILES. Probably from a wrong charge... This needs to be fixed"
                        )
                    temp_G.add_node(
                        str(hash(v[i])),
                        smiles=smiles,
                        #molecule=v[i],
                        E=round((v[i].qm_energies[0]-self.lowest_E)*au2kcal, 1),
                        image=os.path.join(
                            "2D_images", smiles.replace("/", "_") + ".png"
                        ),
                        label=str(round((v[i].qm_energies[0]-self.lowest_E)*au2kcal, 1)),#kcal/mol
                        size= 15,
                        color="green",
                    )
                # else:
                #    temp_G.add_node(
                #        v[i],
                #        # smiles=smiles,
                #        # image=os.path.join(
                #        #    "2D_images", smiles.replace("/", "_") + ".png"
                #        # ),
                #        E=round(float(v[i].qm_energies[0]), 5),
                #        color="red",
                #    )

                if i == 2:
                    temp_G.add_edge(
                        str(hash(v[i])),
                        str(hash(v[i - 2])),
                        #molecule=v[i - 1],
                        label=str(round((v[i-1].qm_energies[0]-self.lowest_E)*au2kcal, 1)),
                    )
        self.G = nx.compose(self.G, temp_G)

        return self.G


def collect(dirs):
    """
    It collects optimized reactants, transition states, products and generates unit reactions.
    """
    result = OrderedDict()
    print("Collecting unit reactions...")
    for d in dirs:
        if d.is_dir() and d.name.isnumeric():
            path_to_rxn = os.path.join(d.path, "gathered/reactions")
            for rxns in os.scandir(path_to_rxn):
                path_to_pathway = os.path.join(rxns.path, "Reaction/pathways/")
                if os.path.exists(path_to_pathway):
                    for unit_rxns in os.scandir(path_to_pathway):
                        if unit_rxns.is_dir() and os.path.exists(
                            os.path.join(unit_rxns, "FS/irc.xyz")
                        ):
                            path_to_FS = os.path.join(unit_rxns, "FS")
                            tar = tarfile.open(
                                os.path.join(path_to_FS, "freezing-string.tar.bz2"),
                                "r:bz2",
                            )
                            tar.extractall(path_to_FS)

                            path_to_react = os.path.join(path_to_FS, "irc_reactant.xyz")
                            path_to_prod = os.path.join(path_to_FS, "irc_product.xyz")
                            path_to_ts = os.path.join(path_to_FS, "irc_transition.xyz")

                            if os.path.exists(path_to_react) and os.path.exists(
                                path_to_prod
                            ):

                                react = Molecule(path_to_react)
                                ts = Molecule(path_to_ts)
                                prod = Molecule(path_to_prod)

                            else:
                                path_to_irc = os.path.join(path_to_FS, "irc.xyz")
                                path_to_ts = os.path.join(path_to_FS, "ts.xyz")
                                react = Molecule(path_to_irc)[0]
                                ts = Molecule(path_to_ts)
                                prod = Molecule(path_to_irc)[-1]

                            react.build_bonds()
                            react.build_topology()
                            ts.build_bonds()
                            ts.build_topology()
                            prod.build_bonds()
                            prod.build_topology()

                            react.qm_energies = [float(react.comms[0].split()[5])]
                            ts.qm_energies = [float(ts.comms[0].split()[-1])]
                            prod.qm_energies = [float(prod.comms[0].split()[5])]

                            key = (
                                d.name
                                + "-"
                                + rxns.name.split("_")[-1]
                                + "_"
                                + unit_rxns.name
                            )
                            result[key] = [react, ts, prod]
    print("Collecting is done. %i unit reactions are ready." % len(result))
    return result

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--charge", type=int, default=0, help="Molecular charges")
    parser.add_argument(
        "--rmsd",
        type=float,
        default=0.05,
        help="If RMSD of two aligned molecule's Cartesian coordinate is better than this threshold, the two molecules are considered different",
    )
    #parser.add_argument(
    #    "--figsize",
    #    type=float,
    #    default=15.0,
    #    help="Size of potential energy surface graphs (inches)",
    #)
    args = parser.parse_args(sys.argv[1:])

    cwd = os.getcwd()
    dirs = os.scandir(cwd)
    rxns = collect(dirs)

    BuildG = BuildGraph(rxns, args.charge)

    G = BuildG.build()

    Energies = list(G.nodes.data('E'))
    Energies.sort(key=lambda a: a[1])
    norm_factor = Energies[-1][-1]

    rgb1 = np.array(colors.to_rgb('#00B2FF'))
    rgb2 = np.array(colors.to_rgb('#FF0000'))

    color_hex = [colors.to_hex((1-a[-1]/norm_factor)*rgb1 + (a[-1]/norm_factor)*rgb2) for a in Energies]

    for i, E in enumerate(Energies):
        G.nodes[E[0]]['color'] = color_hex[i]

    nt = Network('1000px', '1000px')
    nt.from_nx(G)
    nt.toggle_physics(True)
    nt.show('test_color.html')

    print("Done!")


if __name__ == "__main__":
    main()
