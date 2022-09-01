from main import get_conventional_structure, UnitCell
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objs as go
import numpy as np

#from tryout import plotfig



def plot(cell: UnitCell) -> go.Figure:
    atom_number = []

    node_x = []
    node_y = []
    node_z = []

    axis_x = []
    axis_y = []
    axis_z = []
    axes = []

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
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        hovermode="closest",
        height=850
    )

    fig = go.Figure(layout=layout)

    for start, end in cell.edges:
        fig.add_trace(
            go.Scatter3d(
                x=[start[0], end[0], None],
                y=[start[1], end[1], None],
                z=[start[2], end[2], None],
                mode="lines",
                line={
                    "width": 2,
                    "color": "black"
                },
                hoverinfo="none",
                customdata=[1, 1, None]
            )
        )

    origin = np.array([0, 0, 0])
    axes.append((origin, cell.x_axis))
    axes.append((origin, cell.y_axis))
    axes.append((origin, cell.z_axis))
    axes.append((cell.y_axis, cell.y_axis + cell.x_axis))
    # x-axis parallel in z-direction
    axes.append((cell.z_axis, cell.z_axis + cell.x_axis))
    # x-axis parallel in yz-direction
    axes.append((cell.y_axis + cell.z_axis, cell.y_axis + cell.z_axis + cell.x_axis))
    # y-axis parallel in x-direction
    axes.append((cell.x_axis, cell.x_axis + cell.y_axis))
    # y-axis parallel in z-direction
    axes.append((cell.z_axis, cell.z_axis + cell.y_axis))
    # y-axis parallel in xz-direction
    axes.append((cell.x_axis + cell.z_axis, cell.x_axis + cell.z_axis + cell.y_axis))
    # z-axis parallel in x-direction
    axes.append((cell.x_axis, cell.x_axis + cell.z_axis))
    # z-axis parallel in y-direction
    axes.append((cell.y_axis, cell.y_axis + cell.z_axis))
    # z-axis parallel in xy-direction
    axes.append((cell.x_axis + cell.y_axis, cell.x_axis + cell.y_axis + cell.z_axis))

    for axis in axes:
        start, end = axis
        axis_x += [start[0], end[0], None]
        axis_y += [start[1], end[1], None]
        axis_z += [start[2], end[2], None]

    axes_trace = go.Scatter3d(
        x=axis_x,
        y=axis_y,
        z=axis_z,
        mode="lines",
        hoverinfo="none",
        line=dict(color="grey", width=1)
    )

    fig.add_trace(axes_trace)

    for i, coord in enumerate(cell.coords):
        node_x.append(coord[0])
        node_y.append(coord[1])
        node_z.append(coord[2])
        atom_number.append(cell.atoms[i]["number"])

    node_trace = go.Scatter3d(
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

    fig.add_trace(node_trace)

    return fig



dir = "mp-10143"
s = get_conventional_structure(dir)
cell = UnitCell(s)
fig = plot(cell)
#fig = plotfig



app = Dash(__name__)

app.layout = html.Div([
        html.Div(
            dcc.Graph(
                id="graph",
                figure=fig,
                clear_on_unhover=True,
            ),
            id="graph-container",
            style={"border": "2px black solid"},
        ),
        html.Div(
            html.Pre(
                id="helper"
                ),
        ),
        html.Div(
            html.Div([
                dcc.Graph(
                    id="icohp-graph",
                    figure=cell.testgraph,
                ),
                html.Pre(
                    id="bond_strength",
                )
            ]),
            #id="icohp-graph-container",
            id="data-container",
            style={"border":" 2px black solid"},
        ),
    ],
    id="container",
)



@app.callback(
    Output("bond_strength", "children"),
    Output("graph", "figure"),
    Output("icohp-graph", "figure"),
    Input("graph", "hoverData")
)
def edge_hoverevent(hover_data):
    global last_camera_position

    axis = dict(
        showbackground=False,
        showline=False,
        zeroline=False,
        showgrid=False,
        showticklabels=False,
        visible=False,
        title="",
        showspikes=False
    )

    layout = go.Layout(
        showlegend=False,
        scene=dict(
            xaxis=axis,
            yaxis=axis,
        ),
        # margin=dict(t=100),
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        xaxis=dict(visible=False), # TODO: muss am Ende entfernt werden
        yaxis=dict(visible=False), # TODO: muss am Ende entfern werden
        height=250,
        width=250,
    )

    cell.testgraph = go.Figure(layout=layout)
    cell.testgraph.add_trace(go.Scatter(x=[None], y=[None]))

    try:
        data = hover_data["points"][0]["customdata"]
        cell.testgraph.update_traces(x=[1, 2, 3], y=[3, 2, 1])
    except:
        data = "Keine Daten"

    for trace in fig.data:
        if "customdata" in trace:
            trace["line"]["width"] = 2
    if hover_data:
        if "customdata" in hover_data["points"][0]:
            trace_index = hover_data["points"][0]["curveNumber"]
            fig.data[trace_index]["line"]["width"] = 5
            fig.data[trace_index]["opacity"] = 1
    fig.update_layout(scene_camera=last_camera_position)

    return data, fig, cell.testgraph



@app.callback(Output("helper", "children"), Input("graph", "relayoutData"))
def get_current_camera_position(layout_data):
    global last_camera_position
    last_camera_position = dict()
    if layout_data:
        try:
            camera_data = layout_data["scene.camera"]
            up_x = camera_data["up"]["x"]
            up_y = camera_data["up"]["y"]
            up_z = camera_data["up"]["z"]
            center_x = camera_data["center"]["x"]
            center_y = camera_data["center"]["y"]
            center_z = camera_data["center"]["z"]
            eye_x = camera_data["eye"]["x"]
            eye_y = camera_data["eye"]["y"]
            eye_z = camera_data["eye"]["z"]
            last_camera_position = dict(
                up=dict(x=up_x, y=up_y, z=up_z),
                center=dict(x=center_x, y=center_y, z=center_z),
                eye=dict(x=eye_x, y=eye_y, z=eye_z)
            )
        except:
            pass
    return None



if __name__ == '__main__':
    global last_camera_position
    app.run_server(debug=True)