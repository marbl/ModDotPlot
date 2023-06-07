from venv import create
import plotly.express as px
from moddotplot.estimate_identity import (
    binomial_distance,
    poisson_distance,
    containment,
    jaccard,
    get_mods,
    create_coordinates,
)
import itertools
import numpy as np
import dash
from dash import Input, Output, ctx, html, dcc
import dash_daq as daq
from moddotplot.const import COLOR_SCHEMES


def verify_modimizers(sparsity):
    if sparsity % 2 != 0:
        sparsity = sparsity + 1

    sparsity_layers = [sparsity]
    for i in range(4):
        if sparsity_layers[-1] == 1:
            return sparsity_layers
        elif sparsity_layers[-1] % 2 == 1:
            sparsity_layers[-1] = int(sparsity_layers[-1] + 1)
        sparsity_layers.append(int(sparsity_layers[-1] / 2))
    return sparsity_layers


def set_zoom_levels(axis_length, sparsity_layers):
    zoom_levels = {}
    zoom_levels[0] = axis_length
    for i in range(1, len(sparsity_layers)):
        zoom_levels[i] = round(axis_length / pow(2, i))
    return zoom_levels
  

def create_map(kmer_list, resolution, sparsity):
    kmers = get_mods(kmer_list, sparsity)
    kmer_coordinates = create_coordinates2(kmers, resolution)
    return kmer_coordinates


def kmer_coordinates(kmer_list, jacc, resolution, k):
    jaccard_values = []
    for w in itertools.product(kmer_list, repeat=2):
        if jacc:
            j_hat = poisson_distance(containment(w[0], w[1]), k)
            jaccard_values.append(j_hat)
        else:
            j_hat = binomial_distance(containment(w[0], w[1]), k)
            jaccard_values.append(j_hat)

    jaccard_matrix = np.asanyarray(jaccard_values)
    jaccard_matrix.shape = resolution, resolution

    return jaccard_matrix

def kmer_coordinates2(kmer_list1, kmer_list2, jacc, resolution, k):
    jaccard_values = []
    for w in kmer_list1:
        for x in kmer_list2:
            if jacc:
                j_hat = poisson_distance(containment(w, x), k)
                jaccard_values.append(j_hat)
            else:
                j_hat = binomial_distance(containment(w, x), k)
                jaccard_values.append(j_hat)

    jaccard_matrix = np.asanyarray(jaccard_values)
    jaccard_matrix.shape = resolution, resolution

    return jaccard_matrix

def create_coordinates2(kmer_list):
    tmp = []
    for k in kmer_list:
        tmp.append(set())
        for j in k:
            tmp[-1].add(j)

    return tmp


def run_dash(kmer_list1, resolution, sparsity, k, identity, port_number, jacc):
    # Run Dash app
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    current_color = COLOR_SCHEMES["Spectral"]

    # Create initial coordinates & figure
    init = get_mods(kmer_list1, sparsity, resolution)

    new_init = create_coordinates2(init)
    main_level = kmer_coordinates(new_init, False, resolution, k)
    increments = round(len(kmer_list1) / resolution)
    x_ax = [i * increments for i in range(resolution)]
    fig = px.imshow(
        main_level,
        x=x_ax,
        y=x_ax,
        zmin=identity / 100,
        zmax=1,
        height=1000,
        width=1000,
        color_continuous_scale=current_color,
    )

    # Optional layout changes
    fig.update_layout(hovermode="x unified")
    fig.update_xaxes(showspikes=True, spikemode="across")
    fig.update_yaxes(showspikes=True, spikemode="across")
    fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Rockwell")
    )

    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, COLOR_SCHEMES["Spectral"])
    colornames.insert(0, "Spectral")

    # HTML for plot webpage
    app.layout = html.Div(
        [
            daq.ToggleSwitch(
                id="my-toggle-switch",
                value=False,
                color="blue",
                style={
                    "width": "200px",
                    "display": "inline-block",
                    "padding-top": "30px",
                    "padding-left": "80px",
                },
            ),
            html.Div(
                id="my-toggle-switch-output",
                style={
                    "display": "inline-block",
                    "font-family": "Rockwell, Arial, sans-serif",
                    "padding-top": "30px",
                },
            ),
            html.Div(dcc.Graph(id="dotplot", figure=fig)),
        ]
    )

    # Callback activates when panning, zooming,
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Input("dotplot", "relayoutData"),
        Input("my-toggle-switch", "value"),
        Input("dotplot", "figure"),
    )
    def update_x_timeseries(relayoutData, value, figure):
        try:
            if value:
                return figure
            elif relayoutData["xaxis.range[0]"] and not value:
                # TODO: Handle out of bounds error
                if (
                    relayoutData["xaxis.range[0]"] > 0
                    and relayoutData["yaxis.range[1]"] > 0
                ):
                    x_begin = round(relayoutData["xaxis.range[0]"])
                    x_end = round(relayoutData["xaxis.range[1]"])
                    y_begin = round(relayoutData["yaxis.range[0]"])
                    y_end = round(relayoutData["yaxis.range[1]"])
                    amount = x_end - x_begin

                    new_jaccard = []
                    x_mers = kmer_list1[x_begin:x_end]
                    y_mers = kmer_list1[y_end:y_begin]
                    new_map_x = get_mods(x_mers, sparsity, resolution)
                    new_map_y = get_mods(y_mers, sparsity, resolution)

                    x_init = create_coordinates2(new_map_x)
                    y_init = create_coordinates2(new_map_y)
                    updated_jaccard_values = []
                    
                    for w in y_init:
                        for q in x_init:
                            if jacc:
                                updated_jaccard_values.append(
                                    poisson_distance(jaccard(q, w), k)
                                )
                            else:
                                updated_jaccard_values.append(
                                    binomial_distance(containment(q, w), k)
                                )

                    new_jaccard = np.asanyarray(updated_jaccard_values)
                    new_jaccard.shape = resolution, resolution
                    harvest2 = new_jaccard
                    x_ax = [
                        relayoutData["xaxis.range[0]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    y_ax = [
                        relayoutData["yaxis.range[1]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    new_fig = px.imshow(
                        harvest2,
                        x=x_ax,
                        y=y_ax,
                        zmin=identity / 100,
                        zmax=1,
                        height=1000,
                        width=1000,
                        color_continuous_scale=COLOR_SCHEMES["Spectral"],
                    )
                    current_fig = new_fig
                    current_fig.update_layout(hovermode="x unified")
                    current_fig.update_xaxes(showspikes=True, spikemode="across")
                    current_fig.update_yaxes(showspikes=True, spikemode="across")
                    current_fig.update_layout(
                        hoverlabel=dict(
                            bgcolor="white", font_size=16, font_family="Rockwell"
                        )
                    )
                return current_fig
        except TypeError as err:
            # Supress this specific instance of a type error warning, since this will always occur on setup.
            return fig

        except KeyError:
            try:
                # This occurs when the axes are reset. Revert back to the original figure
                if relayoutData["xaxis.autorange"]:
                    return fig
            except KeyError as errw:
                if relayoutData["dragmode"] == "pan":
                    return figure

    # Callback for toggle switch
    @app.callback(
        Output("my-toggle-switch-output", "children"),
        Input("my-toggle-switch", "value"),
    )
    def update_output(value):
        if not value:
            return f"Resolution lock: Off"
        else:
            return f"Resolution lock: On"

    app.run_server(debug=True, use_reloader=False, port=port_number)

def run_dash_pairwise(kmer_list1, kmer_list2, resolution, sparsity, k, identity, port_number, jacc):
    # Run Dash app
    print("using pairwise")
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    current_color = COLOR_SCHEMES["Spectral"]

    # Create initial coordinates & figure
    x_d = get_mods(kmer_list1, sparsity, resolution)
    y_d = get_mods(kmer_list2, sparsity, resolution)

    x_init = create_coordinates2(x_d)
    y_init = create_coordinates2(y_d)

    main_level = kmer_coordinates2(x_init, y_init, False, resolution, k)
    x_increments = round(len(kmer_list1) / resolution)
    y_increments = round(len(kmer_list2) / resolution)
    x_ax = [i * x_increments for i in range(resolution)]
    y_ax = [i * y_increments for i in range(resolution)]
    fig = px.imshow(
        main_level,
        x=x_ax,
        y=y_ax,
        zmin=identity / 100,
        zmax=1,
        height=1000,
        width=1000,
        color_continuous_scale=current_color,
    )

    # Optional layout changes
    fig.update_layout(hovermode="x unified")
    fig.update_xaxes(showspikes=True, spikemode="across")
    fig.update_yaxes(showspikes=True, spikemode="across")
    fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Rockwell")
    )

    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, COLOR_SCHEMES["Spectral"])
    colornames.insert(0, "Spectral")

    # HTML for plot webpage
    app.layout = html.Div(
        [
            daq.ToggleSwitch(
                id="my-toggle-switch",
                value=False,
                color="blue",
                style={
                    "width": "200px",
                    "display": "inline-block",
                    "padding-top": "30px",
                    "padding-left": "80px",
                },
            ),
            html.Div(
                id="my-toggle-switch-output",
                style={
                    "display": "inline-block",
                    "font-family": "Rockwell, Arial, sans-serif",
                    "padding-top": "30px",
                },
            ),
            html.Div(dcc.Graph(id="dotplot", figure=fig)),
        ]
    )

    # Callback activates when panning, zooming,
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Input("dotplot", "relayoutData"),
        Input("my-toggle-switch", "value"),
        Input("dotplot", "figure"),
    )
    def update_x_timeseries(relayoutData, value, figure):
        try:
            if value:
                return figure
            elif relayoutData["xaxis.range[0]"] and not value:
                # TODO: Handle out of bounds error
                if (
                    relayoutData["xaxis.range[0]"] > 0
                    and relayoutData["yaxis.range[1]"] > 0
                ):
                    x_begin = round(relayoutData["xaxis.range[0]"])
                    x_end = round(relayoutData["xaxis.range[1]"])
                    y_begin = round(relayoutData["yaxis.range[0]"])
                    y_end = round(relayoutData["yaxis.range[1]"])
                    amount = x_end - x_begin

                    new_jaccard = []
                    x_mers = kmer_list1[x_begin:x_end]
                    y_mers = kmer_list2[y_end:y_begin]
                    new_map_x = get_mods(x_mers, sparsity, resolution)
                    new_map_y = get_mods(y_mers, sparsity, resolution)

                    x_init = create_coordinates2(new_map_x)
                    y_init = create_coordinates2(new_map_y)
                    updated_jaccard_values = []
                    
                    for w in y_init:
                        for q in x_init:
                            if jacc:
                                updated_jaccard_values.append(
                                    poisson_distance(jaccard(q, w), k)
                                )
                            else:
                                updated_jaccard_values.append(
                                    binomial_distance(containment(q, w), k)
                                )

                    new_jaccard = np.asanyarray(updated_jaccard_values)
                    new_jaccard.shape = resolution, resolution
                    harvest2 = new_jaccard
                    x_ax = [
                        relayoutData["xaxis.range[0]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    y_ax = [
                        relayoutData["yaxis.range[1]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    new_fig = px.imshow(
                        harvest2,
                        x=x_ax,
                        y=y_ax,
                        zmin=identity / 100,
                        zmax=1,
                        height=1000,
                        width=1000,
                        color_continuous_scale=COLOR_SCHEMES["Spectral"],
                    )
                    current_fig = new_fig
                    current_fig.update_layout(hovermode="x unified")
                    current_fig.update_xaxes(showspikes=True, spikemode="across")
                    current_fig.update_yaxes(showspikes=True, spikemode="across")
                    current_fig.update_layout(
                        hoverlabel=dict(
                            bgcolor="white", font_size=16, font_family="Rockwell"
                        )
                    )
                return current_fig
        except TypeError as err:
            # Supress this specific instance of a type error warning, since this will always occur on setup.
            return fig

        except KeyError:
            try:
                # This occurs when the axes are reset. Revert back to the original figure
                if relayoutData["xaxis.autorange"]:
                    return fig
            except KeyError as errw:
                if relayoutData["dragmode"] == "pan":
                    return figure

    # Callback for toggle switch
    @app.callback(
        Output("my-toggle-switch-output", "children"),
        Input("my-toggle-switch", "value"),
    )
    def update_output(value):
        if not value:
            return f"Resolution lock: Off"
        else:
            return f"Resolution lock: On"

    app.run_server(debug=True, use_reloader=False, port=port_number)
