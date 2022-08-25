from pmg import get_conventional_structure, UnitCell
from dash import Dash, dcc, html, Input, Output
import plotly.graph_objs as go
import numpy as np



def generate_plot():
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[1, 2, 3], name="A", line={"width": 1}))
    fig.add_trace(go.Scatter(x=[1, 2, 3], y=[1, 3, 5], name="B", line={"width": 1}))
    return fig



def plot_edges(cell: UnitCell) -> go.Figure:
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

    fig = go.Figure(layout=layout)

    coords = cell.coords.copy()
    for start, end in cell.edges:
        c_a = coords[start]
        c_b = coords[end]
        fig.add_trace(
            go.Scatter3d(
                x=[c_a[0], c_b[0]],
                y=[c_a[1], c_b[1]],
                z=[c_a[2], c_b[2]],
                mode="lines",
                line={
                    "width": 2,
                    "color": "black"
                },
                hoverinfo="none",
            )
        )
    return fig



def plot_axes(cell: UnitCell):
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

    fig = go.Figure(layout=layout)

    axes = []
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

    for start, end in axes:
        fig.add_trace(
            go.Scatter3d(
                x=[start[0], end[0]],
                y=[start[1], end[1]],
                z=[start[2], end[2]],
                mode="lines",
                line={
                    "width": 1,
                    "color": "grey"
                },
                hoverinfo="none",
            )
        )

    return fig



dir = "mp-0"
s = get_conventional_structure(dir)
cell = UnitCell(s)
fig = cell.plot_cell()
#fig_edges = plot_edges(cell)
#fig_axes = plot_axes(cell)



app = Dash(__name__)

app.layout = html.Div([
        html.Div(
            dcc.Graph(
                id="graph",
                figure=fig,
            ),
            id="graph-container",
            style={"border": "2px black solid"},
        ),
        html.Div(
            html.Pre(
                id="helper"
            ),
            id="data-container",
            style={"border":" 2px black solid"},
        ),
        html.Div(
                id="icohp_graph",
        ),
    ],
    id="container",
)

@app.callback(Output("helper", "children"), Input("graph", "hoverData"))
def show_edge_data(hoverData):
    try:
        data = hoverData["points"][0]["customdata"]
        return data
    except:
        #data = None
        return "Keine Daten"
    #return data

@app.callback(Output("icohp_graph", "children"), Input("graph", "hoverData"))
def show_edge_data():
    return cell.testgraph

"""
@app.callback(Output("helper", "children"), Input("edges", "relayoutData"))
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
    else:
        pass
    return None
"""

"""
@app.callback(Output("edges", "figure"), Input("edges", "hoverData"))
def highlight_trace(hover_data):
    global last_camera_position
    # here you set the default settings
    for trace in fig_edges.data:
        trace["line"]["width"] = 1
        trace["opacity"] = 0.5
    if hover_data:
        trace_index = hover_data["points"][0]["curveNumber"]
        fig_edges.data[trace_index]["line"]["width"] = 5
        fig_edges.data[trace_index]["opacity"] = 1
    fig_edges.update_layout(scene_camera=last_camera_position)
    return fig_edges
"""

if __name__ == '__main__':
    #global last_camera_position
    app.run_server(debug=True)