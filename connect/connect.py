#! /usr/bin/env python
import numpy as np
import time
import sys, os, copy
from .molecule import Molecule, TopEqual, MolEqual, Elements
from collections import OrderedDict
import tarfile
import networkx as nx
import matplotlib.pyplot as plt

# import matplotlib.image as mpimg
from xyz2mol import xyz2mol
from rdkit import Chem
from rdkit.Chem import Draw


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
    M = m1 + m2
    M.align()
    rmsd = np.sqrt(np.mean((M[0].xyzs[0] - M[1].xyzs[0]) ** 2))
    return top_equal and rmsd < threshold  # if True else MolEqual(m1, m2)


class BuildGraph(object):
    """
    Use networkx to build a graph
    """

    def __init__(self, rxns, charge):
        self.G = nx.Graph()
        self.rxns = rxns
        self.charge = charge

    def unify_between_lists(self, rxns, ite):
        temp = copy.deepcopy(rxns)
        for i, (k1, v1) in enumerate(rxns.items()):
            for j, (k2, v2) in enumerate(rxns.items()):
                if i < j and i == ite:
                    # temp_list = copy.deepcopy(v2)
                    for k, m1 in enumerate(v1):
                        for l, m2 in enumerate(v2):
                            if l > k and abs(k - l) % 2 == 0 and equal(m1, m2):
                                temp[k2][l] = m1
                                # temp_list[l] = m1
                    # temp[k2] = temp_list
        return temp

    def unify(self, rxns):
        """Pick out different moleucle obejcts with the same connectivity/topology and unify them"""

        temp = copy.deepcopy(rxns)
        ite = len(rxns)
        for i in range(ite - 1):
            temp = self.unify_between_lists(rxns, i)

        unified = copy.deepcopy(temp)
        for k, v in temp.items():
            for i, m1 in enumerate(v):
                for j, m2 in enumerate(v):
                    if i < j and abs(i - j) % 2 == 0 and equal(m1, m2):
                        unified[k][j] = m1

        return unified

    def reduce(self):
        reduced = {}
        for k, v in self.rxns.items():
            reduced_list = []
            for i in range(len(v)):
                if i == 0 or i == len(v) - 1:
                    reduced_list.append(v[i])
                elif i % 3 == 2:
                    E1 = v[i].qm_energies[0]
                    E2 = v[i + 1].qm_energies[0]
                    if E1 < E2:
                        reduced_list.append(v[i])
                    else:
                        reduced_list.append(v[i + 1])
                elif i % 3 == 0:
                    continue
                else:
                    reduced_list.append(v[i])
            reduced[k] = reduced_list
        return reduced

    def build(self):
        temp_G = nx.Graph()
        reduced = self.reduce()
        for k, v in self.unify(reduced).items():
            for i in range(len(v)):
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
                        os.path.join("2D_images", "%s.png" % smiles.replace("/", "_")),
                    )

                else:
                    print(
                        "WARNING: A node is missing SMILES. Probably from a wrong charge... This needs to be fixed"
                    )

                if i % 2 == 1:
                    temp_G.add_node(
                        v[i],
                        smiles=smiles,
                        image=os.path.join(
                            "2D_images", smiles.replace("/", "_") + ".png"
                        ),
                        E=round(float(v[i].qm_energies[0]), 5),
                        color="red",
                    )
                else:
                    temp_G.add_node(
                        v[i],
                        smiles=smiles,
                        image=os.path.join(
                            "2D_images", smiles.replace("/", "_") + ".png"
                        ),
                        E=round(float(v[i].qm_energies[0]), 5),
                        color="green",
                    )
                if i != len(v) - 1:
                    temp_G.add_edge(v[i], v[i + 1])
            self.G = nx.compose(self.G, temp_G)

        # for path in list(nx.connected_components(self.G)):
        #    if len(path) < 6:
        #        for node in path:
        #            self.G.remove_node(node)

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

                            react.qm_energies = [react.comms[0].split()[5]]
                            ts.qm_energies = [ts.comms[0].split()[-1]]
                            prod.qm_energies = [prod.comms[0].split()[5]]

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


def compare_rxns(rxn1, rxn2, threshold):
    """
    Checking to see whether two reaction pathways are the same.

    parameters
    ----------
    rxn1, rxn2: [Molecule objects]
        Molecule objects reactant, TS, and product in a list.

    Return
    ----------
    "True" if they are identical.
    "False" if they are different.
    """
    L1 = len(rxn1)
    L2 = len(rxn2)
    rxns = {L1: rxn1, L2: rxn2}
    # if L1 != L2:
    #     return False
    if L1 == L2:
        equal_count = 0
        for M1, M2 in zip(rxn1, rxn2):
            if equal(M1, M2, threshold):
                equal_count += 1

        if equal_count == L1:
            return True

        equal_count = 0
        for M1, M2 in zip(rxn1[::-1], rxn2):
            if equal(M1, M2, threshold):
                equal_count += 1

        if equal_count == L1:
            return True
    else:
        shorter = rxns[min(rxns)]
        longer = rxns[max(rxns)]
        ite_num = int(abs(L1 / 3 - L2 / 3) + 1)
        equal_count = 0
        for i in range(ite_num):
            for M1, M2 in zip(shorter, longer[i * 3 :]):
                if equal(M1, M2, threshold):
                    equal_count += 1
                    if equal_count == min(L1, L2):
                        return True

            for M1, M2 in zip(shorter, longer[i * 3 :][::-1]):
                if equal(M1, M2, threshold):
                    equal_count += 1
                    if equal_count == min(L1, L2):
                        return True

    return False


def check_repeat(M, threshold):
    """
    Check whether there is a repeating pattern such as BCB in ABCBD. If we allow this pattern, reaction pathway can grow infinitely.
    """
    if len(M) <= 6:
        return False

    for i in range(len(M)):
        if i % 3 == 0:
            for j in range(int(len(M) / 3)):
                if j * 3 > i:
                    if equal(M[i], M[j * 3], threshold) or equal(
                        M[i], M[-1], threshold
                    ):
                        return True
    return False


def filterTS(M_info, E1):
    """ """
    temp = copy.deepcopy(M_info)
    for k, v in M_info.items():
        num = int(len(v) / 3)
        Final_E = float(v[-1].qm_energies[0])
        if Final_E > E1:
            print("Exothermic! We are keeping it.")
            continue
        for i in range(num):
            Rct_E = float(v[i * 3].qm_energies[0])

            TS_E = float(v[i * 3 + 1].qm_energies[0])

            Prd_E = float(v[i * 3 + 2].qm_energies[0])

            if TS_E > Rct_E and TS_E > Prd_E:
                if TS_E > E1:
                    del temp[k]
                    break
            else:
                print("Something is wrong")

    return M_info


def connect_rxns(M_info, iteration=0, outsiders=OrderedDict(), threshold=0.05):
    """
    This function will connect unit reactions.

    Parameters
    ----------
    M_info : OrderedDict
        OrderedDict with frames in key and Molcule objects in a list [reactant, TS, product]

    iteration : int
        Number of iteration

    outsider : OrderedDict
        Reaction pathways that don't make any connections

    Return
    ----------
    final : OrderedDict
        Molecule objects consist of unit reactions
    """

    rxns = OrderedDict()
    temp = copy.deepcopy(M_info)
    connect = 0
    print("-----------------Iteration: %i-----------------" % iteration)
    print("Detecting connection points..")
    for i, (k1, v1) in enumerate(M_info.items()):
        for j, (k2, v2) in enumerate(M_info.items()):
            if j > i:
                if compare_rxns(v1, v2, threshold):
                    # Removing identical reaction pathways
                    if len(v1) >= len(v2) and k2 in temp.keys():
                        del temp[k2]
                    elif len(v1) < len(v2) and k1 in temp.keys():
                        del temp[k1]
                    continue

                reac1 = v1[0]
                prod1 = v1[-1]
                reac2 = v2[0]
                prod2 = v2[-1]
                frm = k1 + "/" + k2
                if equal(reac1, reac2, threshold):
                    M = v2[::-1] + v1
                    if check_repeat(M, threshold):
                        continue
                    else:
                        try:
                            del temp[k1], temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

                elif equal(reac1, prod2, threshold):
                    M = v2 + v1
                    if check_repeat(M, threshold):
                        continue
                    else:
                        try:
                            del temp[k1], temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

                elif equal(prod1, reac2, threshold):
                    M = v1 + v2
                    if check_repeat(M, threshold):
                        continue
                    else:
                        try:
                            del temp[k1], temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

                elif equal(prod1, prod2, threshold):
                    M = v1 + v2[::-1]
                    if check_repeat(M, threshold):
                        continue
                    else:
                        try:
                            del temp[k1], temp[k2]
                        except:
                            pass
                        rxns[frm] = M
                        connect += 1

    outsiders.update(temp)
    print("Number of connections made", connect)
    print("Number of reactions in rxns", len(rxns))
    print("Number of outsiders", len(outsiders))
    if connect == 0:
        rxns.update(outsiders)
    final = copy.deepcopy(rxns)
    print("Numer of reactions saved", len(final))
    print("Cleaning..")
    for i, (k1, v1) in enumerate(rxns.items()):
        for j, (k2, v2) in enumerate(rxns.items()):
            if j > i:
                # if len(v1) == 3 and len(v2) == 3:
                if compare_rxns(v1, v2, threshold):
                    if len(v1) >= len(v2) and k2 in final.keys():
                        del final[k2]
                    elif len(v1) < len(v2) and k1 in final.keys():
                        del final[k1]
    print("After cleaning", len(final))
    iteration += 1
    if connect == 0:
        print("Done! %i reactions detected." % len(final))
        return final
    else:
        return connect_rxns(final, iteration, outsiders, threshold)


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
    parser.add_argument(
        "--figsize",
        type=float,
        default=15.0,
        help="Size of potential energy surface graphs (inches)",
    )
    parser.add_argument(
        "--imgsize",
        type=float,
        default=1.0,
        help="Bigger number will generate smaller 2D molecular images for graphs with more than 3 nodes",
    )
    parser.add_argument(
        "--imgx",
        type=float,
        default=0.0,
        help="Adjust x-axis of 2D molecular images (relative to figsize)",
    )
    parser.add_argument(
        "--imgy",
        type=float,
        default=0.0,
        help="Adjust y-axis of 2D molecular images (relative to figsize)",
    )
    args = parser.parse_args(sys.argv[1:])

    cwd = os.getcwd()
    dirs = os.scandir(cwd)
    result = collect(dirs)
    start = time.time()
    rxns = connect_rxns(result, args.rmsd)

    print("Connecting took %f seconds" % (time.time() - start))

    if not os.path.exists("reactions"):
        os.mkdir("reactions")

    for i, (k, v) in enumerate(rxns.items()):
        Mol = copy.deepcopy(v[0])
        for j in range(len(v) - 1):
            Mol += v[j + 1]
        Mol.write("reactions/%i_%i.xyz" % (i, len(Mol)))

    G = BuildGraph(rxns, args.charge).build()
    subs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    if not os.path.exists("graphs"):
        os.mkdir("graphs")
    for i, G in enumerate(subs):
        pos = nx.spring_layout(G, iterations=1000)
        fig, ax = plt.subplots(figsize=(args.figsize, args.figsize))
        Es = nx.get_node_attributes(G, "E")
        # images = nx.get_node_attributes(G, "image")
        colors = nx.get_node_attributes(G, "color")
        # smiles = nx.get_node_attributes(G, "smiles")
        nx.draw(G, pos=pos, node_color=colors.values(), with_labels=False)
        nodenames = {n: "%f" % Es[n] for n in G.nodes()}
        nx.draw_networkx_labels(G, pos=pos, labels=nodenames, font_weight="bold")
        plt.savefig(os.path.join("graphs", "graph_%i.png" % i))

        tr_figure = ax.transData.transform
        tr_axes = fig.transFigure.inverted().transform
        nnodes = G.number_of_nodes()

        if nnodes == 3:
            img_size = 0.15
        else:
            img_size = 0.15 / (
                np.sqrt(nnodes / 3)) * args.imgsize

        img_center = img_size / 2

        for n in G.nodes:
            if colors[n] != "red":
                xf, yf = tr_figure(pos[n])
                xa, ya = tr_axes((xf, yf))
                a = plt.axes(
                    [
                        xa - img_center + args.imgx,
                        ya - img_center * 2.0 + args.imgy,
                        img_size,
                        img_size,
                    ]
                )
                mol_img = plt.imread(G.nodes[n]["image"])

                imgshape = mol_img.shape
                alpha = np.zeros(imgshape)[:, :, 0]
                alpha[np.where(mol_img[:, :, 0] != 1)] = 1
                mol_img = np.dstack(
                    (mol_img, alpha.reshape(imgshape[0], imgshape[1], 1))
                )
                a.imshow(mol_img)
                a.axis("off")

        plt.savefig("graph_PES_%i.png" % i)

    print("Done!")


if __name__ == "__main__":
    main()
