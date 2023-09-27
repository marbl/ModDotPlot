import plotly.express as px
from moddotplot.estimate_identity import (
    get_mods,
    convert_set,
    convert_set_neighbors,
    pairwise_containment_matrix,
    get_interactive_color,
    get_matching_colors,
    verify_modimizers,
    set_zoom_levels_list,
    make_differences_equal,
    generate_dict_from_list,
    find_value_in_range,
    partition_evenly_spaced_modimizers
)

from moddotplot.static_plots import paired_bed_file

import numpy as np
import dash
from dash import Input, Output, ctx, html, dcc, State
import math
import plotly.graph_objs as go


def run_dash(
    kmer_list1,
    kmer_list2,
    x_name,
    y_name,
    resolution,
    sparsity,
    k,
    identity,
    port_number,
    palette,
    palette_orientation,
    alpha,
    image_pyramid
):
    # Run Dash app
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    app.title = "ModDotPlot"
    # Set color palette
    current_color = get_interactive_color(palette, palette_orientation)
    # Get zooming thresholds, adjust sparsity respectively.
    mod_ranges = verify_modimizers(sparsity, len(image_pyramid))

    mod_thresholds_list = set_zoom_levels_list(round(len(kmer_list1)), mod_ranges)

    important = generate_dict_from_list(mod_thresholds_list)

    image_axes = []
    # Multiply the identity threshold by 100
    for img in image_pyramid:
        img *= 100
        increments = round(len(kmer_list1) / img.shape[0])
        image_axes.append([i * increments for i in range(img.shape[0])])

    main_level = image_pyramid[0]
    main_axis = image_axes[0]

    heatmap = go.Heatmap(z=main_level,
                        zmin=identity,
                        zmax=100,
                        colorscale=current_color,
                        x=main_axis,
                        y=main_axis,
                        showscale=True,  # Shows color scale bar
                        colorbar=dict(title="Identity"),  # Color bar title
                        hoverinfo="x+y+z",  # Hover info to display
                        text=main_level,  # Text to display on hover
                        hovertemplate="X: %{x}<br>Y: %{y}<br><b>Identity: %{z:.2f}</b>",  # Custom hover template
                        name=""
                        )

    fig = go.Figure(data=[heatmap])
    fig.update_xaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True)
    fig.update_xaxes(nticks=10, title_text=x_name, title_font=dict(size=18))

    # Update y-axis properties directly within the heatmap trace
    fig.update_yaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True, autorange='reversed')
    fig.update_yaxes(title_text=y_name, title_font=dict(size=18))

    # Set layout properties
    fig.update_layout(
        height=800,
        width=800,
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica"),
        yaxis_scaleanchor="x",
    )

    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, current_color)
    colornames.insert(0, palette)

    # HTML for plot webpage
    app.layout = html.Div(
        id="main_color",
        style={"display": "flex", "flexDirection": "row", "alignItems": "center", "backgroundColor": "white", "color": "white", "height": "100vh"},
        children=[
            html.Div(
                [
                    dcc.Graph(id="dotplot", figure=fig),
                    html.Div(
                        dcc.RangeSlider(
                            id='threshold-slider',
                            min=identity,
                            max=100,
                            step=0.1,
                            value=[identity, 100], 
                            marks=None,  
                            vertical=True,  
                            verticalHeight=580,  
                            allowCross=False,  
                            pushable=0.1,  
                            tooltip={"placement": "bottom", "always_visible": False},
                            className="custom-slider"
                        ),
                        style={"display": "flex", "justify-content": "center", "padding-top": "110px"},
                    ),
                ],  
                style={"display": "flex", "justify-content": "center"},
            ),
            html.Div(
                style={"width": "330px", "justify-content": "left"},
                children=[
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
                            "padding-top": "140px",
                            "padding-left": "20px",
                            "padding-bot": "130px",
                            "width": "120px",
                            "font-family": "Helvetica, Arial, sans-serif",
                        },  # Added padding to separate the content
                    ),
                    html.Div(
                        html.Button("Toggle Dark Mode", id="dark-mode-toggle"),
                            style={"padding-top": "50px","padding-left": "20px","padding-bot": "30px"}
                    ),
                    html.Div(
                        html.Button("+ Expand Plot", id="zoom-in-button"),
                            style={"padding-top": "10px", "padding-left": "20px","padding-bot": "30px", "width": "160px", "text-align": "left"}
                    ),
                    html.Div(
                        html.Button("- Condense Plot", id="zoom-out-button"),
                            style={"padding-top": "10px", "padding-left": "20px","padding-bot": "30px", "width": "160px", "text-align": "left"}
                    ),
                    html.Div(
                        html.Button("Print plot", id="print"),
                            style={"padding-top": "10px", "padding-left": "20px","padding-bot": "30px", "width": "160px", "text-align": "left"}
                    ),
                    
                    # Placeholder for content
                    html.Div(id="content"),
                    html.Div(
                        "Color Palette",
                        id="color-palette-text",
                        style={
                            "font-family": "Helvetica, Arial, sans-serif",
                            "padding-top": "30px",
                            "padding-left": "20px",  # Added padding to separate the content
                        },
                    ),
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
                    html.Div(
                        f"Alpha: {alpha}",
                            style={"font-family": "Helvetica, Arial, sans-serif", "padding-top": "10px", "padding-left": "20px","padding-bot": "30px", "width": "160px", "text-align": "left"}
                    ),
                ],
            ),
        ],
    )

    #TODO: Install multi tab option for comparative plots
    '''@app.callback(
        Output('dotplot', 'figure', allow_duplicate=True),
        Input('tabs', 'value')
    )
    def render_content(tab):
        if tab == 'tab-1':
            print("one")
            return fig
        elif tab == 'tab-2':
            return fig'''

    #TODO: Install identity slider
    ''' @app.callback(
        Output('dotplot', 'figure'),
        Input('range-slider', 'value'),
        State('dotplot', 'figure')
    )
    def update_slider_output(value, current_figure):
        # Make a copy of the current figure
        #fig = current_figure
        #print(fig)
        # Update the zmin and zmax properties of the figure's layout
        print(fig)
        fig.update_traces(
            zmin=value[0],
            zmax=value[1],
        )

        return fig'''

    @app.callback(
        Output('main_color', 'style', allow_duplicate=True),
        Input('print', 'n_clicks'),
        prevent_initial_call=True,
    )
    def print_test(clickies):
        mod_list = get_mods(kmer_list1, sparsity, resolution)
        windows = partition_evenly_spaced_modimizers(
                    mod_list, len(kmer_list1) + k - 1, resolution
                )
        #print(windows[0][0])
        paired_bed_file(
            windows,
            x_name,
            identity,
            "./",
            False,
            False,
            False,
            "Spectral_11", #TODO: Change this
            "+",
            None,
            None,
            None,
            k,
            False,
            None
        )
        #paired_bed_file()
        dark_mode_style = {
            "backgroundColor": "#161625",
            "color": "white",
            "display": "flex", 
            "height": "100vh", 
        }
        return dark_mode_style

    @app.callback(
        Output('dotplot', 'figure', allow_duplicate=True),
        Input('zoom-in-button', 'n_clicks'),
        Input('zoom-out-button', 'n_clicks'),
        State('dotplot', 'figure'),
        prevent_initial_call=True,
        allow_duplicate=True,
    )
    def zoom_figure(zoom_in_clicks, zoom_out_clicks, current_figure):
        ctx = dash.callback_context

        if ctx.triggered:
            trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

            if trigger_id == "zoom-in-button":
                new_width = current_figure["layout"]['width'] + 100
                new_height = current_figure["layout"]['height'] + 100
            elif trigger_id == "zoom-out-button":
                new_width = max(current_figure["layout"]['width'] - 100, 100)
                new_height = max(current_figure["layout"]['height'] - 100, 100)
            else:
                raise dash.exceptions.PreventUpdate

            fig.update_layout(width=new_width, height=new_height)

        return fig


    # Dark mode
    @app.callback(
        Output("main_color", "style", allow_duplicate=True),
        Input("dark-mode-toggle", "n_clicks"),
        State('dotplot', 'figure'),
    )
    
    def toggle_dark_mode(n_clicks, figure):
        if n_clicks is None:
            n_clicks = 0

        # Toggle dark mode based on the number of clicks
        dark_mode = n_clicks % 2 == 1

        # Define the styles for light and dark modes
        light_mode_style = {
            "backgroundColor": "white",
            "color": "black",
            "display": "flex", 
            "height": "100vh", 
        }

        dark_mode_style = {
            "backgroundColor": "#161625",
            "color": "white",
            "display": "flex", 
            "height": "100vh", 
        }

        '''if dark_mode:
            figure.update_xaxes(showline=True, linewidth=2, linecolor="white", mirror=True)
            figure.update_yaxes(showline=True, linewidth=2, linecolor="white", mirror=True)
        else:
            figure.update_xaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor="black", mirror=True)'''

        # Return the appropriate style based on the current mode
        return dark_mode_style if dark_mode else light_mode_style

    '''@app.callback(
        Output('dotplot', 'figure', allow_duplicate=True),
        [Input('threshold-slider', 'value'),
        Input('gradient-toggle', 'value')]
    )
    def update_heatmap(threshold_range, button_states):
        # Apply the threshold range to the data using conditional indexing
        masked_data = np.where((main_level < threshold_range[0]) | (main_level > threshold_range[1]), 0, main_level)
        print(masked_data[0])
        # Create the heatmap using go.Heatmap
        
        # Adjust the color scale based on the "Keep Original Gradient" button
        if 'keep-original' in button_states:
            zmin = 0
            zmax = 16
        else:
            zmin, zmax = threshold_range
        
        new_heatmap = go.Heatmap(z=masked_data,
                            zmin=zmin,
                            zmax=zmax,
                            x=main_axis,
                            y=main_axis)
        
        # Create a figure with the updated heatmap
        fig = go.Figure(data=[heatmap])
        fig.update_layout(yaxis_scaleanchor="x")
        #print("Is something happening???")
        return fig'''

    # Callback activates when panning, zooming, or changing color
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Output("loading-output-2", "children"),
        Input("dotplot", "relayoutData"),
        Input("dotplot", "figure"),
        Input("color-options", "value"),
        Input("threshold-slider", "value"),
        prevent_initial_call=True,
    )
    def update_dotplot(relayoutData, figure, color, threshold_range):
        current_color = get_interactive_color(get_matching_colors(color), "+")
        new_heatmap = heatmap
        new_heatmap.update(dict(colorscale=current_color))

        masked_data = np.where((new_heatmap["z"] < threshold_range[0]) | (new_heatmap["z"] > threshold_range[1]), 0, new_heatmap["z"])
        fig = go.Figure(data=[new_heatmap])
        fig.update_traces(z=masked_data)
        fig.update_xaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True)
        fig.update_xaxes(nticks=10, title_text=x_name, title_font=dict(size=18))

        # Update y-axis properties directly within the heatmap trace
        fig.update_yaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True, autorange='reversed')
        fig.update_yaxes(title_text=y_name, title_font=dict(size=18))
        # Set layout properties

        # Figure doesn't have height here... how do I fix that? 
        fig.update_layout(
            hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica"),
        )
        
        fig.update_layout(yaxis_scaleanchor="x")
        # Flip y axis since the layout of ModDotPlot it's goes top to bottom
        # TODO: Add option to flip y axis
        if relayoutData is not None:
            if 'xaxis.range[0]' in relayoutData:
                x_start_range = relayoutData.get("xaxis.range[0]")
                y_start_range = relayoutData.get("yaxis.range[1]")
                x_end_range = relayoutData.get("xaxis.range[1]")
                y_end_range = relayoutData.get("yaxis.range[0]")

                # Check that selected range is in bounds, snap to boundary otherwise
                if x_start_range is not None:
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

                    # This function finds the correct level in the image pyramid
                    new_yy = find_value_in_range(amount, important)

                    if new_yy > len(image_pyramid) - 1:
                        x_mers = kmer_list1[x_begin:x_end]
                        y_mers = kmer_list1[y_end:y_begin]

                        x_mods = get_mods(x_mers, sparsity/(2**new_yy) , resolution)
                        y_mods = get_mods(y_mers, sparsity/(2**new_yy), resolution)

                        x_init = convert_set(x_mods)
                        y_init = convert_set(y_mods)
                        x_other = convert_set_neighbors(x_mods, alpha)
                        y_other = convert_set_neighbors(y_mods, alpha)

                        identity_matrix = pairwise_containment_matrix(x_init, y_init, x_other, y_other, identity, k, True)
                        identity_matrix *= 100

                        x_ax = [
                            relayoutData["xaxis.range[0]"] + (i * amount / resolution)
                            for i in range(resolution)
                        ]
                        y_ax = [
                            relayoutData["yaxis.range[1]"] + (i * amount / resolution)
                            for i in range(resolution)
                        ]
                        new_heatmap = go.Heatmap(
                            z=using_matrix,
                            zmin=identity,
                            zmax=100,
                            colorscale=current_color,
                            x=x_ax,
                            y=y_ax,
                            showscale=True,  # Shows color scale bar
                            colorbar=dict(title="Identity"),  # Color bar title
                            hoverinfo="x+y+z",  # Hover info to display
                            text=main_level,  # Text to display on hover
                            hovertemplate="X: %{x}<br>Y: %{y}<br>Identity: %{z:.2f}",
                            name=""
                        )

                        current_fig = go.Figure(data=[new_heatmap])
                        current_fig.update_xaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True)
                        current_fig.update_xaxes(nticks=10, title_text=x_name, title_font=dict(size=18))

                        # Update y-axis properties directly within the heatmap trace
                        current_fig.update_yaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True, autorange='reversed')
                        current_fig.update_yaxes(title_text=y_name, title_font=dict(size=18))
                        current_fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white", font_size=16, font_family="Helvetica"
                            )
                        )
                        current_fig.update_layout(yaxis_scaleanchor="x")
                    
                    else:
                        using_matrix = image_pyramid[new_yy]

                        # Get the coordinates for the better matrix

                        x_start_updated = math.floor(x_begin / len(kmer_list1)*resolution*(2**new_yy))
                        x_end_updated = math.ceil(x_end / len(kmer_list1)*resolution*(2**new_yy))
                        y_start_updated = math.floor(y_end / len(kmer_list1)*resolution*(2**new_yy))
                        y_end_updated = math.ceil(y_begin / len(kmer_list1)*resolution*(2**new_yy))

                        x_equi, y_equi = make_differences_equal(x_start_updated, x_end_updated,y_start_updated,y_end_updated)
                        using_matrix = np.copy(using_matrix[y_start_updated:y_equi, x_start_updated:x_equi])
                        masked_matrix = np.where((using_matrix < threshold_range[0]) | (using_matrix > threshold_range[1]), 0, using_matrix)

                        x_ax = image_axes[new_yy][x_start_updated:x_equi]
                        y_ax = image_axes[new_yy][y_start_updated:y_equi]


                        # Create figure
                        new_heatmap = go.Heatmap(
                            z=masked_matrix,
                            zmin=identity,
                            zmax=100,
                            colorscale=current_color,
                            x=x_ax,
                            y=y_ax,
                            showscale=True,  # Shows color scale bar
                            colorbar=dict(title="Identity"),  # Color bar title
                            hoverinfo="x+y+z",  # Hover info to display
                            text=main_level,  # Text to display on hover
                            hovertemplate="X: %{x}<br>Y: %{y}<br>Identity: %{z:.2f}", 
                            name=""
                        )

                        current_fig = go.Figure(data=[new_heatmap])
                        current_fig.update_xaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True)
                        current_fig.update_xaxes(nticks=10, title_text=x_name, title_font=dict(size=18))

                        # Update y-axis properties directly within the heatmap trace
                        current_fig.update_yaxes(showspikes=True, spikemode="across", ticks="outside", showline=True, linewidth=2, linecolor="black", mirror=True, autorange='reversed')
                        current_fig.update_yaxes(title_text=y_name, title_font=dict(size=18))
                        current_fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white", font_size=16, font_family="Helvetica"
                            )
                        )
                        current_fig.update_layout(yaxis_scaleanchor="x")
                        return current_fig, ""
                    
            elif 'xaxis.autorange' in relayoutData:
                # Triggers when double-clicked or reset axes
                return fig, ""

            #TODO: disable Autosizing
            elif 'autosize' in relayoutData:
                # Autosize should be disabled: if not it gets handled here
                return figure, ""

            else:
                # Happens when user selects pan mode. 
                return figure, ""
        else:
            #figure.update_layout(coloraxis=dict(colorscale=current_color))
            # This occurs at startup and when selecting color. Return the base figure with color
            return fig, ""

        
    app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            /* Custom CSS for the RangeSlider */
            .custom-slider .rc-slider-track {
                background-color: gray;  /* Change the highlight color to gray */
            }

            .custom-slider .rc-slider-dot {
                display: none;  /* Hide the dots */
            }

            .custom-slider .rc-slider-handle {
                border: 2px solid black;  /* Change the handle border color to gray */
                font-family: Helvetica, Arial, sans-serif;
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''
    
    app.run_server(debug=True, use_reloader=False, port=port_number)