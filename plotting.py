import os
import plotly.graph_objs as go

from pymatgen.analysis.graphs import Structure, StructureGraph
from structuregraph_helpers.create import get_structure_graph
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer



def plot(sg: StructureGraph) -> None:
    node_x = []
    node_y = []
    node_z = []

    edge_x = []
    edge_y = []
    edge_z = []

    atom_number = []

    coords = sg.structure.frac_coords

    for i, _ in enumerate(coords):
        c = coords[i]

        node_x.append(c[0])
        node_y.append(c[1])
        node_z.append(c[2])

        atom_number.append(sg.structure[i].specie.number)

    for start, end, data in sg.graph.edges(data=True):
        start_c = coords[start]
        end_c = coords[end] + data["to_jimage"]

        edge_x += [start_c[0], end_c[0], None]
        edge_y += [start_c[1], end_c[1], None]
        edge_z += [start_c[2], end_c[2], None]

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

    fig = go.Figure(data=data, layout=layout)
    fig.show()



#dir = "mp-1779/"
dir = "mp-0/"
ext = os.path.expanduser("~/automationresults/")
filepath = ext + dir + "POSCAR"

struct = Structure.from_file(filepath)
struct = SpacegroupAnalyzer(struct)
struct = struct.get_conventional_standard_structure()

scaling = [[2,0,0], [0,1,0], [0,0,1]]
print(struct)
struct.make_supercell(scaling_matrix=scaling)
print(struct)