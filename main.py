import os
import random

import plotly.graph_objs as go
import plotly.express as px
import numpy as np
import networkx as nx

from structuregraph_helpers.create import get_structure_graph
from plotly.graph_objs import Figure
from pymatgen.analysis.graphs import Structure, StructureGraph
from structuregraph_helpers.create import get_structure_graph
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from random import randrange
from itertools import permutations, product, combinations_with_replacement
from pymatgen.analysis.local_env import (
    BrunnerNN_relative,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    MinimumDistanceNN,
    NearNeighbors,
    VoronoiNN,
)



class UnitCell:

    def __init__(self, structure: Structure):
        self.testgraph = px.line(x=[None], y=[None])
        self.structuregraph = get_structure_graph(structure, method="minimaldistance")
        self.number_individual_atoms = len(structure.frac_coords)
        self.a = structure.lattice.a
        self.b = structure.lattice.b
        self.c = structure.lattice.c
        self.alpha = structure.lattice.alpha
        self.beta = structure.lattice.beta
        self.gamma = structure.lattice.gamma

        self.x_axis = np.array(
            [
                self.a,
                0,
                0
            ]
        )
        self.y_axis = np.array(
            [
                self.b * np.cos(np.deg2rad(self.gamma)),
                self.b * np.sin(np.deg2rad(self.gamma)),
                0
            ]
        )
        self.z_axis = np.array(
            [
                self.c * np.cos(np.deg2rad(self.beta)),
                self.c * np.cos(np.deg2rad(self.alpha)),
                self.c * np.sin(np.deg2rad(self.gamma))
            ]
        )

        self.frac_coords = list(structure.frac_coords.copy())
        self.cart_crystal_axis_mat = np.stack((self.x_axis, self.y_axis, self.z_axis), axis=-1)
        self.coords = []
        self.atoms = dict()
        self.eq_atoms = dict()
        new_coords = []
        for i, coord in enumerate(structure.frac_coords.copy()):
            self.coords.append(np.dot(self.cart_crystal_axis_mat, coord))
            self.atoms[i] = dict()
            self.atoms[i]["element"] = structure[i].specie
            self.atoms[i]["number"] = structure[i].specie.number
            self.eq_atoms[i] = []
            indices0 = []
            indices1 = []
            for j, c in enumerate(coord):
                if c < 0.01:
                    indices0.append(j)
                elif c > 0.99:
                    indices1.append(j)
            n0 = len(indices0)
            n1 = len(indices1)
            add0 = [np.zeros(shape=3)]
            add1 = [np.zeros(shape=3)]
            for j in range(1, n0 + 1):
                v = [1] * j + [0] * (n0 - n1 - j)
                ps0 = set(permutations(v))
                for p0 in ps0:
                    a0 = np.zeros(shape=3)
                    for index, permutation in zip(indices0, p0):
                        a0[index] = permutation
                    add0.append(a0)
            for j in range(1, n1 + 1):
                v = [-1] * j + [0] * (n1 - n0 - j)
                ps1 = set(permutations(v))
                for p1 in ps1:
                    a1 = np.zeros(shape=3)
                    for index, permutation in zip(indices1, p1):
                        a1[index] = permutation
                    add1.append(a1)
            adds = product(add0, add1)
            for a0, a1 in adds:
                if np.linalg.norm(a0 + a1) > 0:
                    new_frac_coord = coord.copy() + a0 + a1
                    self.frac_coords.append(new_frac_coord)
                    new_coord = np.array(coord.copy()) + a0 + a1
                    new_coord = np.dot(self.cart_crystal_axis_mat, new_coord)
                    new_coords.append(new_coord)
                    m = self.number_individual_atoms - 1 + len(new_coords)
                    self.atoms[m] = dict()
                    self.atoms[m]["element"] = structure[i].specie
                    self.atoms[m]["number"] = structure[i].specie.number
                    self.eq_atoms[i].append(m)
        self.coords += new_coords

        self.edges = []
        limit = 0.01
        for node1, node2, data in self.structuregraph.graph.edges(data=True):
            start = self.frac_coords[node1]
            end = self.frac_coords[node2] + data["to_jimage"]
            start_coord = np.dot(self.cart_crystal_axis_mat, start)
            end_coord = np.dot(self.cart_crystal_axis_mat, end)
            d0 = np.linalg.norm(start_coord - end_coord)
            if (-limit <= end[0] <= limit+1) and \
               (-limit <= end[1] <= limit+1) and \
               (-limit <= end[2] <= limit+1):
                self.edges.append((start_coord, end_coord))
            for eq_atom1 in self.eq_atoms[node1]:
                shift = self.frac_coords[eq_atom1] - self.frac_coords[node1]
                new_end = end + shift
                if (-limit <= new_end[0] <= limit+1) and \
                   (-limit <= new_end[1] <= limit+1) and \
                   (-limit <= new_end[2] <= limit+1):
                    start_coord = self.coords[eq_atom1]
                    end_coord = np.dot(self.cart_crystal_axis_mat, new_end)
                    d1 = np.linalg.norm(start_coord - end_coord)
                    if abs(d0 - d1) <= limit:
                        self.edges.append((start_coord, end_coord))
            # TODO: Müssen hier noch die anderen jeweils äquivalenten Atome durchgegangen werden, oder reicht das schon?


    def plot_cell(self):
        atom_number = []

        node_x = []
        node_y = []
        node_z = []

        axis_x = []
        axis_y = []
        axis_z = []
        axes = []

        edge_x = []
        edge_y = []
        edge_z = []

        origin = np.array([0, 0, 0])
        axes.append((origin, self.x_axis))
        axes.append((origin, self.y_axis))
        axes.append((origin, self.z_axis))
        axes.append((self.y_axis, self.y_axis + self.x_axis))
        # x-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.x_axis))
        # x-axis parallel in yz-direction
        axes.append((self.y_axis + self.z_axis, self.y_axis + self.z_axis + self.x_axis))
        # y-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.y_axis))
        # y-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.y_axis))
        # y-axis parallel in xz-direction
        axes.append((self.x_axis + self.z_axis, self.x_axis + self.z_axis + self.y_axis))
        # z-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.z_axis))
        # z-axis parallel in y-direction
        axes.append((self.y_axis, self.y_axis + self.z_axis))
        # z-axis parallel in xy-direction
        axes.append((self.x_axis + self.y_axis, self.x_axis + self.y_axis + self.z_axis))
        
        for axis in axes:
            start, end = axis
            axis_x += [start[0], end[0], None]
            axis_y += [start[1], end[1], None]
            axis_z += [start[2], end[2], None]

        for i, coord in enumerate(self.coords):
            node_x.append(coord[0])
            node_y.append(coord[1])
            node_z.append(coord[2])
            atom_number.append(self.atoms[i]["number"])

        #for start, end in self.edges:
        #    start_coord = self.coords[start]
        #    end_coord = self.coords[end]
        #    edge_x += [start_coord[0], end_coord[0], None]
        #    edge_y += [start_coord[1], end_coord[1], None]
        #    edge_z += [start_coord[2], end_coord[2], None]

        for start, end in self.edges:
            edge_x += [start[0], end[0], None]
            edge_y += [start[1], end[1], None]
            edge_z += [start[2], end[2], None]


        trace0 = go.Scatter3d(
            x=axis_x,
            y=axis_y,
            z=axis_z,
            mode="lines",
            hoverinfo="none",
            line=dict(color="grey", width=1)
        )

        trace1 = go.Scatter3d(
            x=node_x,
            y=node_y,
            z=node_z,
            mode="markers",
            hoverinfo="none",
            marker=dict(
                symbol="circle",
                size=6,
                color=atom_number,
                colorscale="Viridis",
                line=dict(color="rgb(50,50,50)", width=0.5),
            ),
        )

        trace2 = go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            customdata=None,
            hoverinfo="none",
            line=dict(color="black", width=2),
        )
        
        axis = dict(
            showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title="",
            showspikes=False
        )

        layout = go.Layout(
            showlegend=False,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
            ),
            #margin=dict(t=100),
            margin=dict(
                l=20,
                r=20,
                b=10,
                t=10,
            ),
            hovermode="closest",
            height=850
        )

        data = []
        data.append(trace0)
        data.append(trace1)
        data.append(trace2)

        fig = go.Figure(data=data, layout=layout)
        #fig.show()
        #return fig


    def plot_extended_cell(self, add_cells = None) -> None:
        if len(add_cells) == 0 or add_cells is None:
            self.plot_cell()
            return

        atom_number = []

        node_x = []
        node_y = []
        node_z = []

        axis_x = []
        axis_y = []
        axis_z = []
        axes = []

        edge_x = []
        edge_y = []
        edge_z = []

        origin = np.array([0, 0, 0])
        axes.append((origin, self.x_axis))
        axes.append((origin, self.y_axis))
        axes.append((origin, self.z_axis))
        axes.append((self.y_axis, self.y_axis + self.x_axis))
        # x-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.x_axis))
        # x-axis parallel in yz-direction
        axes.append((self.y_axis + self.z_axis, self.y_axis + self.z_axis + self.x_axis))
        # y-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.y_axis))
        # y-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.y_axis))
        # y-axis parallel in xz-direction
        axes.append((self.x_axis + self.z_axis, self.x_axis + self.z_axis + self.y_axis))
        # z-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.z_axis))
        # z-axis parallel in y-direction
        axes.append((self.y_axis, self.y_axis + self.z_axis))
        # z-axis parallel in xy-direction
        axes.append((self.x_axis + self.y_axis, self.x_axis + self.y_axis + self.z_axis))

        for axis in axes:
            start, end = axis
            axis_x += [start[0], end[0], None]
            axis_y += [start[1], end[1], None]
            axis_z += [start[2], end[2], None]
            for add_cell in add_cells:
                new_start = start.copy()
                new_end = end.copy()
                n = add_cell[0] + 1
                for i in range(1, n):
                    a = add_cell[1]
                    b = add_cell[2]
                    c = add_cell[3]
                    new_start = new_start + a * self.x_axis + b * self.y_axis + c * self.z_axis
                    new_end = new_end + a * self.x_axis + b * self.y_axis + c * self.z_axis
                    axis_x += [new_start[0], new_end[0], None]
                    axis_y += [new_start[1], new_end[1], None]
                    axis_z += [new_start[2], new_end[2], None]

        for i, coord in enumerate(self.coords):
            node_x.append(coord[0])
            node_y.append(coord[1])
            node_z.append(coord[2])
            atom_number.append(self.atoms[i]["number"])
            for add_cell in add_cells:
                new_coord = coord.copy()
                n = add_cell[0] + 1
                for j in range(1, n):
                    a = add_cell[1]
                    b = add_cell[2]
                    c = add_cell[3]
                    new_coord = new_coord + a * self.x_axis + b * self.y_axis + c * self.z_axis
                    node_x.append(new_coord[0])
                    node_y.append(new_coord[1])
                    node_z.append(new_coord[2])
                    atom_number.append(self.atoms[i]["number"])

        for start, end in self.edges:
            edge_x += [start[0], end[0], None]
            edge_y += [start[1], end[1], None]
            edge_z += [start[2], end[2], None]

        trace0 = go.Scatter3d(
            x=axis_x,
            y=axis_y,
            z=axis_z,
            mode="lines",
            hoverinfo="none",
            line=dict(color="lightgrey", width=1)
        )

        trace1 = go.Scatter3d(
            x=node_x,
            y=node_y,
            z=node_z,
            mode="markers",
            hoverinfo="none",
            marker=dict(
                symbol="circle",
                size=6,
                color=atom_number,
                colorscale="Viridis",
                line=dict(color="rgb(50,50,50)", width=0.5),
            ),
        )

        trace2 = go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode="lines",
            hoverinfo="none",
            line=dict(color="black", width=2),
        )

        axis = dict(
            showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title="",
            showspikes=False
        )

        layout = go.Layout(
            showlegend=False,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
            ),
            margin=dict(t=100),
            hovermode="closest",
        )

        data = []
        data.append(trace0)
        data.append(trace1)
        data.append(trace2)

        fig = go.Figure(data=data, layout=layout)
        fig.show()


    def plot_supercell(self, scale: int = None)-> None:
        if scale == None:
            self.plot_cell()
            return

        atom_number = []

        node_x = []
        node_y = []
        node_z = []

        axis_x = []
        axis_y = []
        axis_z = []
        axes = []

        origin = np.array([0, 0, 0])
        axes.append((origin, self.x_axis))
        axes.append((origin, self.y_axis))
        axes.append((origin, self.z_axis))
        axes.append((self.y_axis, self.y_axis + self.x_axis))
        # x-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.x_axis))
        # x-axis parallel in yz-direction
        axes.append((self.y_axis + self.z_axis, self.y_axis + self.z_axis + self.x_axis))
        # y-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.y_axis))
        # y-axis parallel in z-direction
        axes.append((self.z_axis, self.z_axis + self.y_axis))
        # y-axis parallel in xz-direction
        axes.append((self.x_axis + self.z_axis, self.x_axis + self.z_axis + self.y_axis))
        # z-axis parallel in x-direction
        axes.append((self.x_axis, self.x_axis + self.z_axis))
        # z-axis parallel in y-direction
        axes.append((self.y_axis, self.y_axis + self.z_axis))
        # z-axis parallel in xy-direction
        axes.append((self.x_axis + self.y_axis, self.x_axis + self.y_axis + self.z_axis))

        combs = combinations_with_replacement([-1, 0, 1], 3)
        perms = []
        for comb in combs:
            ps = set(permutations(comb))
            for p in ps:
                perms.append(p)

        for axis in axes:
            start, end = axis
            axis_x += [start[0], end[0], None]
            axis_y += [start[1], end[1], None]
            axis_z += [start[2], end[2], None]
            for x, y, z in perms:
                new_axis = axis + (scale-1) * (x * self.x_axis + y * self.y_axis + z * self.z_axis)
                start, end = new_axis
                axis_x += [start[0], end[0], None]
                axis_y += [start[1], end[1], None]
                axis_z += [start[2], end[2], None]

        for i, coord in enumerate(self.coords):
            node_x.append(coord[0])
            node_y.append(coord[1])
            node_z.append(coord[2])
            atom_number.append(self.atoms[i]["number"])
            for x, y, z in perms:
                new_coord = coord + (scale-1) * (x * self.x_axis + y * self.y_axis + z * self.z_axis)
                node_x.append(new_coord[0])
                node_y.append(new_coord[1])
                node_z.append(new_coord[2])
                atom_number.append(self.atoms[i]["number"])

        trace0 = go.Scatter3d(
            x=axis_x,
            y=axis_y,
            z=axis_z,
            mode="lines",
            hoverinfo="none",
            line=dict(color="lightgrey", width=1)
        )

        trace1 = go.Scatter3d(
            x=node_x,
            y=node_y,
            z=node_z,
            mode="markers",
            hoverinfo="none",
            marker=dict(
                symbol="circle",
                size=6,
                color=atom_number,
                colorscale="Viridis",
                line=dict(color="rgb(50,50,50)", width=0.5),
            ),
        )

        axis = dict(
            showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title="",
            showspikes=False
        )

        layout = go.Layout(
            showlegend=False,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
            ),
            margin=dict(t=100),
            hovermode="closest",
        )

        data = []
        data.append(trace0)
        data.append(trace1)

        fig = go.Figure(data=data, layout=layout)
        fig.show()



def plot_structure_graph(sg: StructureGraph) -> Figure:
    node_x = []
    node_y = []
    node_z = []

    edge_x = []
    edge_y = []
    edge_z = []

    axis_x = []
    axis_y = []
    axis_z = []

    atom_number = []

    # setup for axes of unit cell
    a = sg.structure.lattice.a
    b = sg.structure.lattice.b
    c = sg.structure.lattice.c
    alpha = sg.structure.lattice.alpha
    beta = sg.structure.lattice.beta
    gamma = sg.structure.lattice.gamma
    alpha_rad = np.deg2rad(alpha)
    beta_rad = np.deg2rad(beta)
    gamma_rad = np.deg2rad(gamma)
    axes = []
    origin = np.array([0, 0, 0])

    # x-axis
    x = np.array([a, 0, 0])
    axes.append((origin, x))
    # y-axis
    y = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
    axes.append((origin, y))
    # z-axis
    z = np.array([c * np.cos(beta_rad), c * np.cos(alpha_rad), c * np.sin(gamma_rad)])
    axes.append((origin, z))
    # x-axis parallel in y-direction
    axes.append((y, y+x))
    # x-axis parallel in z-direction
    axes.append((z, z+x))
    # x-axis parallel in yz-direction
    axes.append((y+z, y+z+x))
    # y-axis parallel in x-direction
    axes.append((x, x+y))
    # y-axis parallel in z-direction
    axes.append((z, z+y))
    # y-axis parallel in xz-direction
    axes.append((x+z, x+z+y))
    # z-axis parallel in x-direction
    axes.append((x, x+z))
    # z-axis parallel in y-direction
    axes.append((y, y+z))
    # z-axis parallel in xy-direction
    axes.append((x+y, x+y+z))

    for axis in axes:
        start, end = axis
        axis_x += [start[0], end[0], None]
        axis_y += [start[1], end[1], None]
        axis_z += [start[2], end[2], None]

    coords = sg.structure.frac_coords
    axis_trans_mat = np.stack((x, y, z), axis=-1)
    edge_coords = []
    for i, coord in enumerate(coords):
        new_coord = np.dot(axis_trans_mat, coord)
        node_x.append(new_coord[0])
        node_y.append(new_coord[1])
        node_z.append(new_coord[2])
        atom_number.append(sg.structure[i].specie.number)
        #############################
        edge_coords.append(new_coord)
        #############################

    for i, coord in enumerate(coords):
        indices0 = []
        indices1 = []
        for j, c in enumerate(coord):
            if c < 0.01:
                indices0.append(j)
            elif c > 0.99:
                indices1.append(j)
        n0 = len(indices0)
        n1 = len(indices1)
        add0 = [np.zeros(shape=3)]
        add1 = [np.zeros(shape=3)]
        """
        # same functionality as below
        for n, factor, indices, add in ((n0, 1, indices0, add0), (n1, -1, indices1, add1)):
            for j in range(1, n + 1):
                v = [factor] * j + [0] * (n0 - n1 - j)
                ps = set(permutations(v))
                for p in ps:
                    a = np.zeros(shape=3)
                    for index, permutation in zip(indices, p):
                        a[index] = permutation
                    add.append(a)
        """
        for j in range(1, n0 + 1):
            v = [1] * j + [0] * (n0 - n1 - j)
            ps0 = set(permutations(v))
            for p0 in ps0:
                a0 = np.zeros(shape=3)
                for index, permutation in zip(indices0, p0):
                    a0[index] = permutation
                add0.append(a0)
        for j in range(1, n1 + 1):
            v = [-1] * j + [0] * (n1 - n0 - j)
            ps1 = set(permutations(v))
            for p1 in ps1:
                a1 = np.zeros(shape=3)
                for index, permutation in zip(indices1, p1):
                    a1[index] = permutation
                add1.append(a1)
        adds = product(add0, add1)
        for a0, a1 in adds:
            if np.linalg.norm(a0 + a1) > 0:
                new_coord = np.array(coord.copy()) + a0 + a1
                new_coord = np.dot(axis_trans_mat, new_coord)
                node_x.append(new_coord[0])
                node_y.append(new_coord[1])
                node_z.append(new_coord[2])
                atom_number.append(sg.structure[i].specie.number)

    for start, end, data in sg.graph.edges(data=True):
        a, b, c = data["to_jimage"]
        if a == 0 and b == 0 and c == 0:
            start_c = edge_coords[start]
            end_c = edge_coords[end] + data["to_jimage"]
            edge_x += [start_c[0], end_c[0], None]
            edge_y += [start_c[1], end_c[1], None]
            edge_z += [start_c[2], end_c[2], None]

    trace0 = go.Scatter3d(
        x=axis_x,
        y=axis_y,
        z=axis_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="lightgrey", width=1)
    )


    trace1 = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="black", width=2),
    )


    trace2 = go.Scatter3d(
        x=node_x,
        y=node_y,
        z=node_z,
        mode="markers",
        hoverinfo="none",
        marker=dict(
            symbol="circle",
            size=6,
            color=atom_number,
            colorscale="Viridis",
            line=dict(color="rgb(50,50,50)", width=0.5),
        ),
    )

    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        title="",
        showspikes=False
    )

    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=dict(axis),
            yaxis=dict(axis),
            zaxis=dict(axis),
        ),
        margin=dict(t=100),
        hovermode="closest",
    )

    data = []
    data.append(trace2)
    data.append(trace1)
    data.append(trace0)

    fig = go.Figure(data=data, layout=layout)
    return fig



def get_random_crystalsystem_example(crystalsystem: str) -> Structure:
    filepath = os.path.join(os.path.expanduser("~/RAW_files_phonon/"), crystalsystem)
    mps = os.listdir(filepath)
    rand_index = randrange(len(mps))
    mp = mps[rand_index]
    filepath = os.path.join(filepath, mp, "POSCAR")
    print(f"{crystalsystem.upper()} EXAMPLE: {mp}")
    return Structure.from_file(filepath)


def get_structure(dir: str, crystalsystem: str = None) -> Structure:
    if crystalsystem:
        filepath = os.path.join(os.path.expanduser("~/RAW_files_phonon/"), crystalsystem, dir, "POSCAR")
    else:
        filepath = os.path.join(os.path.expanduser("~/automationresults"), dir, "POSCAR")
    s = Structure.from_file(filepath)
    return s


def get_conventional_structure(dir: str, crystalsystem: str = None) -> Structure:
    if crystalsystem:
        filepath = os.path.join(os.path.expanduser("~/RAW_files_phonon/"), crystalsystem, dir, "POSCAR")
    else:
        filepath = os.path.join(os.path.expanduser("~/automationresults"), dir, "POSCAR")
    s = Structure.from_file(filepath)
    sga = SpacegroupAnalyzer(s)
    cs = sga.get_conventional_standard_structure()
    return cs



#crystalsystem = "cubic"
#s = get_random_crystalsystem_example(crystalsystem)
#sga = SpacegroupAnalyzer(s)
#s = sga.get_conventional_standard_structure()

#dir = "mp-553374"
#dir = "mp-997086"
#dir = "mp-684690"
#dir = "mp-7773"
#dir = "mp-6450"

#dir = "mp-567841"
#dir = "mp-29398"
#s = get_conventional_structure(dir, crystalsystem)

dir = "mp-10143/"
#dir = "mp-0/"
#s = get_conventional_structure(dir)
s = get_structure(dir)

#sg = get_structure_graph(s, method="minimaldistance")

#cell = UnitCell(s)
#cell.plot_cell()
#vecs = [
#    [1, 1, 0, 0],
#    [1, -1, 0, 0],
#    [1, 0, 1, 0],
#    [1, 0, -1, 0],
#    [1, 0, 0, 1],
#    [1, 0, 0, -1]
#]
#cell.plot_extended_cell(vecs)
#cell.plot_supercell(2)