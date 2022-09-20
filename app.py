import os

import plotly.graph_objs as go

from helper import get_structure_plot, get_dummy_cohp_plot, get_chosen_structure
from dash import Dash, dcc, html, Input, Output



dir = "mp-10143/"
#dir = "mp-510401"
#dir = "mp-2384"
path = get_chosen_structure(dir)
#path = get_random_structure()
path_to_poscar = os.path.join(path, "POSCAR")
path_to_charge = os.path.join(path, "CHARGE.lobster")
path_to_icobilist = os.path.join(path, "ICOBILIST.lobster")
path_to_icooplist = os.path.join(path, "ICOOPLIST.lobster")
path_to_icohplist = os.path.join(path, "ICOHPLIST.lobster")
path_to_cohpcar = os.path.join(path, "COHPCAR.lobster")
path_to_madelung = os.path.join(path, "MadelungEnergies.lobster")

fig = get_structure_plot(
    path_to_poscar,
    path_to_charge,
    path_to_icobilist,
    path_to_icooplist,
    path_to_icohplist,
    path_to_cohpcar,
    path_to_madelung
)

cohp_plot = get_dummy_cohp_plot()


app = Dash(__name__)

app.layout = html.Div([
        html.Div(
            html.Div( # dieser Hilfscontainer wird nur benötigt, wenn alle Container sichtbare runde Ränder haben
                dcc.Graph(
                    id="structuregraph",
                    figure=fig,
                    clear_on_unhover=True,
                ),
                id="helper-container",
            ), # siehe Kommentar über diesem
            className="container-class",
            id="structuregraph-container",
        ),
        html.Div(
            html.Pre(
                id="helper"
            ),
        ),
        html.Div(
            dcc.Graph(
                id="plot",
                figure=cohp_plot,
            ),
            className="container-class",
            id="plot-container",
        ),
        html.Div(
            html.Table([
                html.Tr([
                    html.Td(
                        "Bond Length",
                        className="property-name"
                    ),
                    html.Td(
                        id="bond_length",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOBI",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOBI",
                        className = "property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOOP",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOOP",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOHP",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOHP",
                        className="property-data"
                    )
                ]),
                html.Tr([
                    html.Td(
                        "ICOHP Bonding Perc",
                        className="property-name"
                    ),
                    html.Td(
                        id="ICOHP_bonding_perc",
                        className="property-data"
                    )
                ]),
            ]),
            className="container-class",
            id="data-container",
        ),
    ],
    id="container",
)



@app.callback(
    Output("bond_length", "children"),
    Output("ICOBI", "children"),
    Output("ICOOP", "children"),
    Output("ICOHP", "children"),
    Output("ICOHP_bonding_perc", "children"),
    Output("structuregraph", "figure"),
    Output("plot", "figure"),
    Input("structuregraph", "hoverData")
)
def edge_hoverevent(hover_data):
    global last_camera_position

    """
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
    """

    layout = go.Layout(
        showlegend=False,
        margin=dict(
            l=20,
            r=20,
            b=10,
            t=10,
        ),
        height=440,
        width=600,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False)
    )

    cohp = go.Figure(layout=layout)

    try:
        cohp_data, bond_length, icobi, icoop, icohp, icohp_bonding_perc = hover_data["points"][0]["customdata"]
        bond_length = f"{bond_length} {chr(8491)}"
        #icobi = ohne Einheit
        #icoop = Anteil Elektronen??
        icohp = f"{icohp} eV"
        icohp_bonding_perc = f"{icohp_bonding_perc*100} %"

        cohp.add_trace(
            go.Scatter(
                x=cohp_data[0],
                y=cohp_data[1],
                line=dict(color="red")
            )
        )
        cohp.add_trace(
            go.Scatter(
                x=[min(cohp_data[0]), max(cohp_data[0])],
                y=[0, 0],
                mode="lines",
                line=dict(
                    width=1,
                    color="black",
                    dash="dash"
                )
            )
        )

        cohp.update_xaxes(
            visible=True,
            title="COHP",
            showline=True,
            linewidth=1.5,
            linecolor="black",
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
        )
        cohp.update_yaxes(
            visible=True,
            title="E - Efermi [eV]",
            showline=True,
            linewidth=1.5,
            linecolor="black",
        )
    except:
        cohp.add_trace(go.Scatter(x=[None], y=[None], line=dict(color="red")))
        bond_length = "-"
        icobi = "-"
        icoop = "-"
        icohp = "-"
        icohp_bonding_perc = "-"

    for trace in fig.data:
        if "customdata" in trace:
            trace["line"]["width"] = 2
    if hover_data:
        if "customdata" in hover_data["points"][0]:
            trace_index = hover_data["points"][0]["curveNumber"]
            fig.data[trace_index]["line"]["width"] = 5
            fig.data[trace_index]["opacity"] = 1
    fig.update_layout(scene_camera=last_camera_position)

    return bond_length, icobi, icoop, icohp, icohp_bonding_perc, fig, cohp



@app.callback(Output("helper", "children"), Input("structuregraph", "relayoutData"))
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