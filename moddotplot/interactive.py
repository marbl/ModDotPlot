from nis import match
from venv import create
import plotly.express as px
from moddotplot.estimate_identity import (
    get_mods,
    convert_set,
    convert_set_neighbors,
    initial_identity_matrix,
    updated_identity_matrix,
    get_interactive_color,
    get_matching_colors,
    verify_modimizers,
    set_zoom_levels,
)
import numpy as np
import dash
from dash import Input, Output, ctx, html, dcc
import dash_daq as daq


def run_dash(
    kmer_list1,
    resolution,
    sparsity,
    k,
    identity,
    port_number,
    jacc,
    palette,
    palette_orientation,
    alpha,
):
    # Run Dash app
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    # Set color palette
    current_color = get_interactive_color(palette, palette_orientation)

    # Get zooming thresholds, adjust sparsity respectively.
    mod_ranges = verify_modimizers(sparsity)
    mod_thresholds = set_zoom_levels(round(len(kmer_list1)), mod_ranges)
    alpha = alpha

    # Create initial coordinates & figure
    mod_list = get_mods(kmer_list1, sparsity, resolution)
    mod_set = convert_set(mod_list)
    main_level = initial_identity_matrix(mod_set, resolution, k)
    main_level *= 100
    increments = round(len(kmer_list1) / resolution)
    x_ax = [i * increments for i in range(resolution)]

    # Create figure
    fig = px.imshow(
        main_level,
        x=x_ax,
        y=x_ax,
        zmin=identity,
        zmax=100,
        height=1000,
        width=1000,
        color_continuous_scale=current_color,
        labels=dict(color="Identity"),
    )

    # Cross layout on hovering dotplot
    fig.update_layout(hovermode="x unified")
    fig.update_xaxes(showspikes=True, spikemode="across")
    fig.update_yaxes(showspikes=True, spikemode="across")
    fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica")
    )

    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, current_color)
    colornames.insert(0, palette)

    # HTML for plot webpage
    app.layout = html.Div(
        style={"display": "flex", "height": "100vh"},
        children=[
            # Dotplot div
            html.Div(
                dcc.Graph(id="dotplot", figure=fig),
                style={"flex": "1", "overflow": "auto"},
            ),
            html.Div(
                style={"flex": "1", "overflow": "auto"},
                children=[
                    # Div for resolution lock switch
                    daq.ToggleSwitch(
                        id="my-toggle-switch",
                        value=False,
                        color="blue",
                        style={
                            "width": "80px",
                            "padding-top": "120px",
                            "padding-left": "20px",
                        },
                    ),
                    # Div for resolution lock text
                    html.Div(
                        id="my-toggle-switch-output",
                        style={
                            "font-family": "Helvetica, Arial, sans-serif",
                            "padding-top": "10px",
                            "padding-left": "20px",
                        },
                    ),
                    # Div for loading screen upon refresh
                    html.Div(
                        [
                            dcc.Loading(
                                id="loading-2",
                                children=[html.Div(id="loading-output-2")],
                                type="circle",
                            )
                        ],
                        style={
                            "padding-top": "50px",
                            "padding-left": "20px",
                            "padding-bot": "30px",
                            "width": "120px",
                            "font-family": "Helvetica, Arial, sans-serif",
                        },  # Added padding to separate the content
                    ),
                    html.Div(
                        "Color Palette",
                        id="color-palette-text",
                        style={
                            "font-family": "Helvetica, Arial, sans-serif",
                            "padding-top": "30px",
                            "padding-left": "20px",  # Added padding to separate the content
                        },
                    ),
                    # TODO: Download plot button
                    # html.Button('Download as png', id='btn-nclicks-1', n_clicks=0, style={ 'justify-content': 'center', 'width': 100, 'height': 100}),
                    # html.Div(id='container-button-timestamp'),
                    html.Div(
                        style={"display": "flex", "height": "100vh"},
                        children=[
                            html.Div(
                                [
                                    dcc.Dropdown(
                                        [
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/Spectral_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Spectral",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Spectral",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Sequential Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Sequential",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Blues_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Blues",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Blues",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/BuGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BuGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BuGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/BuPu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BuPu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BuPu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/GnBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "GnBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "GnBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Greens_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Greens",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Green",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Greys_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Greys",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Greys",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/OrRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "OrRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "OrRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Oranges_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Oranges",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Oranges",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuBuGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuBuGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuBuGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Purples_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Purples",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Purples",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/RdPu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdPu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdPu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Reds_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Reds",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Reds",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlGnBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlGnBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlGnBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlOrBr_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlOrBr",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlOrBr",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlOrRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlOrRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlOrRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Divergent Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Divergent",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/BrBG_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BrBG",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BrBG",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PRGn_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PRGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PRGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PiYG_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PiYG",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PiYG",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PuOr_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuOr",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuOr",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdBu_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdGy_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdGy",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdGy",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdYlBu_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdYlBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdYlBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdYlGn_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdYlGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdYlGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Qualitative Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Qualitative",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Accent_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Accent",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Accent",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Dark2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Dark2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Dark2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Paired_12_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Paired",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Paired",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Pastel1_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Pastel1",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Pastel1",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Pastel2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Pastel2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Pastel2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set1_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set1",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set1",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set3_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set3",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set3",
                                            },
                                        ],
                                        value=palette.split("_")[0],
                                        id="color-options",
                                    ),
                                    html.Div(id="output-div"),
                                ],
                                style={
                                    "padding-top": "10px",
                                    "padding-left": "20px",
                                    "width": "250px",
                                    "font-family": "Helvetica, Arial, sans-serif",
                                },  # Added padding to separate the content
                            )
                        ],
                    ),
                ],
            ),
        ],
    )

    # TODO: download plot code
    '''@app.callback(
        Output("loading-output-2", "children", allow_duplicate=True),
        Input('btn-nclicks-1', 'n_clicks'),
        Input('dotplot', 'figure'),
    )
    def displayClick(btn1,figure):
        if "btn-nclicks-1" == ctx.triggered_id:
            windows = partition_evenly_spaced_modimizers(
                        mod_list, len(kmer_list1) + k - 1, resolution
                    )
            # use as input for pllot download
            print(figure["data"][0]["z"])
    
        paired_bed_file(
                    windows,
                    seq_list[i],
                    args.identity,
                    args.sparsity,
                    args.output,
                    len(k_list[i]) + args.kmer - 1,
                    args.no_bed,
                    args.no_plot,
                    args.no_hist,
                    args.palette,
                    args.palette_orientation,
                    args.width,
                    args.height,
                    args.dpi,
                    args.kmer,
                    args.bin_freq,
                )
        return ""'''

    # Callback activates when panning or zooming
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Output("loading-output-2", "children"),
        Input("dotplot", "relayoutData"),
        Input("my-toggle-switch", "value"),
        Input("dotplot", "figure"),
        Input("color-options", "value"),
    )
    def update_x_timeseries(relayoutData, value, figure, color):
        current_color = get_interactive_color(get_matching_colors(color), "+")
        fig.update_layout(coloraxis=dict(colorscale=current_color))
        try:
            if value:
                return figure, ""
            elif relayoutData["xaxis.range[0]"] and not value:
                # Flip y axis since it's going top to bottom
                x_start_range = relayoutData.get("xaxis.range[0]")
                y_start_range = relayoutData.get("yaxis.range[1]")
                x_end_range = relayoutData.get("xaxis.range[1]")
                y_end_range = relayoutData.get("yaxis.range[0]")

                # Check that selected range is in bounds, snap to boundary otherwise
                if x_start_range is not None and not value:
                    if x_start_range < 0:
                        print("x axis out of bounds! Shifting to x=0")
                        relayoutData["xaxis.range[0]"] = 0
                        x_start_range = 0
                    if y_start_range < 0:
                        print("y axis out of bounds! Shifting to y=0")
                        relayoutData["yaxis.range[1]"] = 0
                        y_start_range = 0
                    if x_end_range > len(kmer_list1):
                        print(f"x axis out of bounds! Shifting to x={len(kmer_list1)}")
                        relayoutData["xaxis.range[1]"] = len(kmer_list1)
                        x_end_range = len(kmer_list1)
                    if y_end_range > len(kmer_list1):
                        print(f"y axis out of bounds! Shifting to y={len(kmer_list1)}")
                        relayoutData["yaxis.range[0]"] = len(kmer_list1)
                        y_end_range = len(kmer_list1)

                if x_start_range >= 0 and y_start_range >= 0:
                    x_begin = round(relayoutData["xaxis.range[0]"])
                    x_end = round(relayoutData["xaxis.range[1]"])

                    y_begin = round(relayoutData["yaxis.range[0]"])
                    y_end = round(relayoutData["yaxis.range[1]"])

                    amount = x_end - x_begin

                    new_sparsity = 1
                    for element in mod_thresholds:
                        if mod_thresholds[element] > amount:
                            new_sparsity = element

                    # Switch order of y
                    x_mers = kmer_list1[x_begin:x_end]
                    y_mers = kmer_list1[y_end:y_begin]
                    x_mods = get_mods(x_mers, new_sparsity, resolution)
                    y_mods = get_mods(y_mers, new_sparsity, resolution)

                    x_init = convert_set(x_mods)
                    y_init = convert_set(y_mods)
                    x_other = convert_set_neighbors(x_mods, alpha)
                    y_other = convert_set_neighbors(y_mods, alpha)

                    identity_matrix = updated_identity_matrix(
                        x_init, y_init, x_other, y_other, resolution, k
                    )
                    identity_matrix *= 100

                    x_ax = [
                        relayoutData["xaxis.range[0]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    y_ax = [
                        relayoutData["yaxis.range[1]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    new_fig = px.imshow(
                        identity_matrix,
                        x=x_ax,
                        y=y_ax,
                        zmin=identity,
                        zmax=100,
                        height=1000,
                        width=1000,
                        color_continuous_scale=current_color,
                    )
                    current_fig = new_fig
                    current_fig.update_layout(hovermode="x unified")
                    current_fig.update_xaxes(showspikes=True, spikemode="across")
                    current_fig.update_yaxes(showspikes=True, spikemode="across")
                    current_fig.update_layout(
                        hoverlabel=dict(
                            bgcolor="white", font_size=16, font_family="Helvetica"
                        )
                    )
                return current_fig, ""
        except TypeError as err:
            # Supress this specific instance of a type error warning, since this will always occur on setup.
            return fig, ""

        except KeyError:
            try:
                # This occurs when the axes are reset. Revert back to the original figure
                if relayoutData["xaxis.autorange"]:
                    return fig, ""
            except KeyError as errw:
                if relayoutData["dragmode"] == "pan":
                    figure.update_layout(coloraxis=dict(colorscale=current_color))
                    return figure, ""

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


def run_dash_pairwise(
    kmer_list1,
    kmer_list2,
    resolution,
    sparsity,
    k,
    identity,
    port_number,
    jacc,
    palette,
    palette_orientation,
    alpha,
):
    # Run Dash app
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    current_color = get_interactive_color(palette, palette_orientation)

    # Create initial coordinates & figure
    x_d = get_mods(kmer_list1, sparsity, resolution)
    y_d = get_mods(kmer_list2, sparsity, resolution)

    x_init = convert_set(x_d)
    y_init = convert_set(y_d)
    x_neighbor = convert_set_neighbors(x_d, alpha)
    y_neighbor = convert_set_neighbors(y_d, alpha)

    # main_level = kmer_coordinates2(x_init, y_init, False, resolution, k)
    main_level = updated_identity_matrix(
        x_init, y_init, x_neighbor, y_neighbor, resolution, k
    )
    main_level *= 100
    x_increments = round(len(kmer_list1) / resolution)
    y_increments = round(len(kmer_list2) / resolution)
    x_ax = [i * x_increments for i in range(resolution)]
    y_ax = [i * y_increments for i in range(resolution)]
    fig = px.imshow(
        main_level,
        x=x_ax,
        y=y_ax,
        zmin=identity,
        zmax=100,
        height=1000,
        width=1000,
        color_continuous_scale=current_color,
        labels=dict(color="Identity"),
    )

    # Optional layout changes
    fig.update_layout(hovermode="x unified")
    fig.update_xaxes(showspikes=True, spikemode="across")
    fig.update_yaxes(showspikes=True, spikemode="across")
    fig.update_layout(
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica")
    )

    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, current_color)
    colornames.insert(0, palette)

    # HTML for plot webpage
    app.layout = html.Div(
        style={"display": "flex", "height": "100vh"},
        children=[
            # Dotplot div
            html.Div(
                dcc.Graph(id="dotplot", figure=fig),
                style={"flex": "1", "overflow": "auto"},
            ),
            html.Div(
                style={"flex": "1", "overflow": "auto"},
                children=[
                    # Div for resolution lock switch
                    daq.ToggleSwitch(
                        id="my-toggle-switch",
                        value=False,
                        color="blue",
                        style={
                            "width": "80px",
                            "padding-top": "120px",
                            "padding-left": "20px",
                        },
                    ),
                    # Div for resolution lock text
                    html.Div(
                        id="my-toggle-switch-output",
                        style={
                            "font-family": "Helvetica, Arial, sans-serif",
                            "padding-top": "10px",
                            "padding-left": "20px",
                        },
                    ),
                    # Div for loading screen upon refresh
                    html.Div(
                        [
                            dcc.Loading(
                                id="loading-2",
                                children=[html.Div(id="loading-output-2")],
                                type="circle",
                            )
                        ],
                        style={
                            "padding-top": "50px",
                            "padding-left": "20px",
                            "padding-bot": "30px",
                            "width": "120px",
                            "font-family": "Helvetica, Arial, sans-serif",
                        },  # Added padding to separate the content
                    ),
                    html.Div(
                        "Color Palette",
                        id="color-palette-text",
                        style={
                            "font-family": "Helvetica, Arial, sans-serif",
                            "padding-top": "30px",
                            "padding-left": "20px",  # Added padding to separate the content
                        },
                    ),
                    # TODO: Download plot button
                    # html.Button('Download as png', id='btn-nclicks-1', n_clicks=0, style={ 'justify-content': 'center', 'width': 100, 'height': 100}),
                    # html.Div(id='container-button-timestamp'),
                    html.Div(
                        style={"display": "flex", "height": "100vh"},
                        children=[
                            html.Div(
                                [
                                    dcc.Dropdown(
                                        [
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/Spectral_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Spectral",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Spectral",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Sequential Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Sequential",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Blues_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Blues",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Blues",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/BuGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BuGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BuGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/BuPu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BuPu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BuPu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/GnBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "GnBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "GnBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Greens_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Greens",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Green",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Greys_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Greys",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Greys",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/OrRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "OrRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "OrRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Oranges_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Oranges",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Oranges",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuBuGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuBuGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuBuGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/PuRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Purples_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Purples",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Purples",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/RdPu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdPu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdPu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/Reds_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Reds",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Reds",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlGn_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlGnBu_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlGnBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlGnBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlOrBr_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlOrBr",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlOrBr",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/sequential/img/YlOrRd_9_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "YlOrRd",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "YlOrRd",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Divergent Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Divergent",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/BrBG_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "BrBG",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "BrBG",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PRGn_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PRGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PRGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PiYG_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PiYG",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PiYG",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/PuOr_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "PuOr",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "PuOr",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdBu_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdGy_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdGy",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdGy",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdYlBu_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdYlBu",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdYlBu",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/diverging/img/RdYlGn_11_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "RdYlGn",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "RdYlGn",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Span(
                                                            "Qualitative Color Palettes",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Qualitative",
                                                "disabled": True,
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Accent_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Accent",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Accent",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Dark2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Dark2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Dark2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Paired_12_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Paired",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Paired",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Pastel1_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Pastel1",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Pastel1",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Pastel2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Pastel2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Pastel2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set1_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set1",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set1",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set2_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set2",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set2",
                                            },
                                            {
                                                "label": html.Span(
                                                    [
                                                        html.Img(
                                                            src="https://jiffyclub.github.io/palettable/colorbrewer/qualitative/img/Set3_8_discrete.png",
                                                            height=10,
                                                            width=100,
                                                        ),
                                                        html.Span(
                                                            "Set3",
                                                            style={
                                                                "font-size": 15,
                                                                "padding-left": 10,
                                                            },
                                                        ),
                                                    ],
                                                    style={
                                                        "align-items": "center",
                                                        "justify-content": "center",
                                                    },
                                                ),
                                                "value": "Set3",
                                            },
                                        ],
                                        value=palette.split("_")[0],
                                        id="color-options",
                                    ),
                                    html.Div(id="output-div"),
                                ],
                                style={
                                    "padding-top": "10px",
                                    "padding-left": "20px",
                                    "width": "250px",
                                    "font-family": "Helvetica, Arial, sans-serif",
                                },  # Added padding to separate the content
                            )
                        ],
                    ),
                ],
            ),
        ],
    )

    # Callback activates when panning, zooming,
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Output("loading-output-2", "children"),
        Input("dotplot", "relayoutData"),
        Input("my-toggle-switch", "value"),
        Input("dotplot", "figure"),
        Input("color-options", "value"),
    )
    def update_x_timeseries(relayoutData, value, figure, color):
        current_color = get_interactive_color(get_matching_colors(color), "+")
        fig.update_layout(coloraxis=dict(colorscale=current_color))
        try:
            if value:
                return figure, ""
            elif relayoutData["xaxis.range[0]"] and not value:
                # Flip y axis since it's going top to bottom
                x_start_range = relayoutData.get("xaxis.range[0]")
                y_start_range = relayoutData.get("yaxis.range[1]")
                x_end_range = relayoutData.get("xaxis.range[1]")
                y_end_range = relayoutData.get("yaxis.range[0]")

                # Check that selected range is in bounds, snap to boundary otherwise
                if x_start_range is not None and not value:
                    if x_start_range < 0:
                        print("x axis out of bounds! Shifting to x=0")
                        relayoutData["xaxis.range[0]"] = 0
                        x_start_range = 0
                    if y_start_range < 0:
                        print("y axis out of bounds! Shifting to y=0")
                        relayoutData["yaxis.range[1]"] = 0
                        y_start_range = 0
                    if x_end_range > len(kmer_list1):
                        print(f"x axis out of bounds! Shifting to x={len(kmer_list1)}")
                        relayoutData["xaxis.range[1]"] = len(kmer_list1)
                        x_end_range = len(kmer_list1)
                    if y_end_range > len(kmer_list1):
                        print(f"y axis out of bounds! Shifting to y={len(kmer_list1)}")
                        relayoutData["yaxis.range[0]"] = len(kmer_list1)
                        y_end_range = len(kmer_list1)

                if x_start_range >= 0 and y_start_range >= 0:
                    x_begin = round(relayoutData["xaxis.range[0]"])
                    x_end = round(relayoutData["xaxis.range[1]"])

                    y_begin = round(relayoutData["yaxis.range[0]"])
                    y_end = round(relayoutData["yaxis.range[1]"])

                    amount = x_end - x_begin
                    y_mers = kmer_list2[x_begin:x_end]
                    x_mers = kmer_list1[y_end:y_begin]
                    x_mods = get_mods(x_mers, sparsity, resolution)
                    y_mods = get_mods(y_mers, sparsity, resolution)

                    x_init = convert_set(x_mods)
                    y_init = convert_set(y_mods)
                    x_other = convert_set_neighbors(x_mods, alpha)
                    y_other = convert_set_neighbors(y_mods, alpha)

                    identity_matrix = updated_identity_matrix(
                        x_init, y_init, x_other, y_other, resolution, k
                    )
                    identity_matrix *= 100

                    x_ax = [
                        relayoutData["xaxis.range[0]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    y_ax = [
                        relayoutData["yaxis.range[1]"] + (i * amount / resolution)
                        for i in range(resolution)
                    ]
                    new_fig = px.imshow(
                        identity_matrix,
                        x=x_ax,
                        y=y_ax,
                        zmin=identity,
                        zmax=100,
                        height=1000,
                        width=1000,
                        color_continuous_scale=current_color,
                    )
                    current_fig = new_fig
                    current_fig.update_layout(hovermode="x unified")
                    current_fig.update_xaxes(showspikes=True, spikemode="across")
                    current_fig.update_yaxes(showspikes=True, spikemode="across")
                    current_fig.update_layout(
                        hoverlabel=dict(
                            bgcolor="white", font_size=16, font_family="Helvetica"
                        )
                    )
                return current_fig, ""
        except TypeError as err:
            # Supress this specific instance of a type error warning, since this will always occur on setup.
            return fig, ""

        except KeyError:
            try:
                # This occurs when the axes are reset. Revert back to the original figure
                if relayoutData["xaxis.autorange"]:
                    return fig, ""
            except KeyError as errw:
                if relayoutData["dragmode"] == "pan":
                    return figure, ""

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
