import plotly.express as px
from moddotplot.estimate_identity import (
    getInteractiveColor,
    getMatchingColors,
    makeDifferencesEqual,
    generateDictionaryFromList,
    findValueInRange,
    convertMatrixToBed,
)

import numpy as np
import dash
from dash import Input, Output, html, dcc, State
import math
import plotly.graph_objs as go
import logging
import os

# Prevent HTTP protocol requests from showing up in terminal
log = logging.getLogger("werkzeug")
log.setLevel(logging.ERROR)


def find_closest_elements(value, sorted_list):
    # Initialize variables to store the indices of the closest elements
    closest_index1 = None
    closest_index2 = None

    # Initialize a variable to store the minimum difference
    min_difference = float("inf")

    # Iterate through the sorted list
    for i in range(len(sorted_list)):
        # Calculate the absolute difference between the current element and the target value
        current_difference = abs(sorted_list[i] - value)

        # Update indices and minimum difference if a closer element is found
        if current_difference < min_difference:
            min_difference = current_difference
            closest_index1 = i

            # Check if there is a previous element to avoid index out of range
            if i > 0:
                closest_index2 = i - 1
    if closest_index1 == None:
        return closest_index2, closest_index2
    elif closest_index2 == None:
        return closest_index1, closest_index1
    if closest_index1 > closest_index2:
        return closest_index2, closest_index1
    else:
        return closest_index1, closest_index2


def run_dash(matrices, metadata, axes, sparsity, identity, port_number, output_dir):
    # Run Dash app
    app = dash.Dash(__name__, prevent_initial_callbacks="initial_duplicate")
    app.title = "ModDotPlot"
    print(
        f"{app.title} interactive mode is successfully running on http://127.0.0.1:{port_number}/ \n"
    )
    # Initialize first sequence
    image_pyramid = matrices[0]
    current_metadata = metadata[0]
    x_name = metadata[0]["x_name"]
    y_name = metadata[0]["y_name"]

    # Initialize color palette
    palette = "Spectral_11"
    palette_orientation = "+"
    current_color = getInteractiveColor(palette, palette_orientation)

    # Initialize axes and titles for all images
    titles = []
    for i in range(len(metadata)):
        titles.append(metadata[i]["title"])

    # Get zooming thresholds, adjust sparsity respectively.
    def halving_sequence(size, start):
        sequence = [start]
        for _ in range(1, size):
            start /= 2
            sequence.append(start)
        return sequence

    # print(current_metadata)
    mod_thresholds_list = halving_sequence(
        len(current_metadata["sparsities"]), current_metadata["x_size"]
    )
    # print(mod_thresholds_list)

    # print(current_metadata["min_window_size"]* current_metadata["resolution"])
    # print(current_metadata["max_window_size"])
    numo = round(
        math.log2(
            current_metadata["max_window_size"] / current_metadata["min_window_size"]
        )
        + 1
    )
    # print(numo)
    important = generateDictionaryFromList(mod_thresholds_list)
    # print(f"this is imprtant: {important}")

    main_level = image_pyramid[0]
    main_x_axis = axes[0][0]
    main_y_axis = axes[0][1]
    main_x_axis_np = np.array(main_x_axis)

    # TODO: modify value here
    main_x_axis_np += 3000

    # Modify text so that hover shows interval format
    # TODO: Add intervals in a way that's speedy
    '''X, Y = np.meshgrid(main_x_axis, main_y_axis)
    X = np.delete(X, -1, axis=1)
    X, Y = np.meshgrid(main_x_axis, main_y_axis)
    X = np.delete(X, -1, axis=1)
    X_str = np.array([[f"{str(val)}-{str(val + current_metadata['max_window_size'])}" for val in row] for row in X])
    Y = np.delete(Y, -1, axis=1)
    Y_str = np.array([[f"{str(val)}-{str(val + current_metadata['max_window_size'])}" for val in row] for row in Y])
    hover_template_text = "%{text}<br>%{customdata}<br>Identity: %{z:.1f}"'''
    hover_template_text = "Identity: %{z:.2f}"
    # Create initial heatmap
    heatmap = go.Heatmap(
        z=main_level,
        zmin=identity,
        zmax=100,
        colorscale=current_color,
        x=main_x_axis,
        y=main_y_axis,
        showscale=True,
        colorbar=dict(title="Identity"),
        hoverinfo="all",
        hovertemplate=hover_template_text,
        name="",
        x0=0,
        dx=current_metadata["max_window_size"],
        xtype="scaled",
        y0=0,
        dy=current_metadata["max_window_size"],
        ytype="scaled",
    )

    fig = go.Figure(data=[heatmap])
    fig.update_xaxes(
        showspikes=True,
        spikemode="across",
        ticks="outside",
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    fig.update_xaxes(
        nticks=10, title_text=metadata[0]["x_name"], title_font=dict(size=18)
    )

    # Update y-axis properties directly within the heatmap trace
    fig.update_yaxes(
        showspikes=True,
        spikemode="across",
        ticks="outside",
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    fig.update_yaxes(title_text=metadata[0]["y_name"], title_font=dict(size=18))
    fig_title = ""

    # TODO: Fine tuning on the names
    title_size = 28
    if metadata[0]["self"]:
        fig_title = f"Self-Identity Plot: {x_name}"
        if len(x_name) > 22:
            title_size = 20
    else:
        fig_title = f"Comparative Plot: {x_name} vs. {y_name}"
        if len(x_name) + len(y_name) > 22:
            title_size = 18

    # Set layout properties
    fig.update_layout(
        height=800,
        width=800,
        hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica"),
        yaxis_scaleanchor="x",
        title=fig_title,
        title_font=dict(size=title_size, family="Helvetica, Arial, sans-serif"),
        title_x=0.5,
        title_y=0.95,
    )
    colorscales = px.colors.named_colorscales()
    colornames = px.colors.named_colorscales()

    colorscales.insert(0, current_color)
    colornames.insert(0, palette)

    # HTML for plot webpage
    app.layout = html.Div(
        id="main_color",
        style={
            "alignItems": "center",
            "backgroundColor": "white",
            "color": "black",
            "height": "100vh",
        },
        children=[
            html.Link(rel="icon", href="/assets/favicon.png"),
            html.Div(
                [
                    dcc.Graph(id="dotplot", figure=fig),
                    html.Div(
                        dcc.RangeSlider(
                            id="threshold-slider",
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
                            className="custom-slider",
                        ),
                        style={
                            "display": "flex",
                            "justifyContent": "center",
                            "paddingTop": "110px",
                        },
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    dcc.Loading(
                                        id="loading-2",
                                        children=[html.Div(id="loading-output-2")],
                                        type="circle",
                                    )
                                ],
                                style={
                                    "paddingBottom": "40px",
                                    "paddingLeft": "40px",
                                    "width": "100px",
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                },  # Added padding to separate the content
                            ),
                            html.Div(
                                id="invisible-div",
                                children=0,
                                style={"display": "none"},
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        "Currently shown plot:",
                                        style={
                                            "marginRight": "10px",
                                            "padding-left": "20px",
                                        },
                                    ),
                                    dcc.Dropdown(
                                        id="matrix-dropdown",
                                        options=[
                                            {"label": f"{title}", "value": f"{title}"}
                                            for title in titles  # Iterate over each title
                                        ],
                                        value=(
                                            titles[0] if titles else None
                                        ),  # Set default value based on the length of matrices
                                        clearable=False,  # Prevent dropdown from clearing values,
                                        style={"width": "300px"},
                                    ),
                                    html.Div(id="output-container"),
                                    dcc.Store(id="matrices-store", data=len(matrices)),
                                ],
                                style={
                                    "display": "none" if len(titles) < 2 else "block",
                                    "width": "fit-content"
                                    * 2,  # Set width to fit content
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                    "paddingTop": "10px",
                                    "paddingLeft": "25px",
                                    "paddingBottom": "30px",
                                    "textAlign": "left",
                                },
                            ),
                            html.Div(
                                html.Button(
                                    "Save Matrix to Bed File",
                                    id="save-bed",
                                    n_clicks=0,
                                    disabled=False,
                                    style={
                                        "marginLeft": "45px",
                                    },
                                ),
                            ),
                            html.Div(
                                children=[
                                    f"Current Window Size: {current_metadata['max_window_size']}"
                                ],
                                id="window-div",
                                style={
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                    "paddingTop": "10px",
                                    "paddingLeft": "45px",
                                    "paddingBot": "30px",
                                    "width": "220px",
                                    "textAlign": "left",
                                },
                            ),
                            html.Div(
                                f"Minimum Window Size: {current_metadata['min_window_size']}",
                                style={
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                    "paddingTop": "10px",
                                    "paddingLeft": "45px",
                                    "paddingBot": "30px",
                                    "width": "220px",
                                    "textAlign": "left",
                                },
                            ),
                            dcc.Checklist(
                                id="gradient-toggle",
                                options=[
                                    {
                                        "label": "  Readjust Gradient",
                                        "value": "keep-original",
                                    }
                                ],
                                value=[],
                                style={
                                    "paddingTop": "10px",
                                    "paddingLeft": "24px",
                                    "paddingBottom": "30px",
                                    "width": "260px",
                                    "textAlign": "left",
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                },
                            ),
                            html.Div(id="content"),
                            html.Div(
                                "Color Palette",
                                id="color-palette-text",
                                style={
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                    "paddingLeft": "45px",  # Added padding to separate the content
                                },
                            ),
                            html.Div(
                                style={"display": "flex", "height": "10vh"},
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                                                        "fontSize": 15,
                                                                        "paddingLeft": 10,
                                                                    },
                                                                ),
                                                            ],
                                                            style={
                                                                "alignItems": "center",
                                                                "justifyContent": "center",
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
                                            "paddingTop": "10px",
                                            "paddingLeft": "20px",
                                            "width": "250px",
                                            "fontFamily": "Helvetica, Arial, sans-serif",
                                        },  # Added padding to separate the content
                                    )
                                ],
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        "Coordinate Log:",
                                        id="coordinate-text",
                                        style={
                                            "fontFamily": "Helvetica, Arial, sans-serif",
                                            "paddingLeft": "45px",  # Added padding to separate the content
                                            "paddingBottom": "10px",
                                        },
                                    ),
                                    dcc.Loading(
                                        id="loading",
                                        type="default",
                                        children=[
                                            dcc.Textarea(
                                                id="text-input",
                                                value="",
                                                style={
                                                    "padding": "10px",
                                                    "height": "40px",
                                                    "width": "62%",
                                                    "margin-left": "22px",
                                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                                },
                                            ),
                                            html.Button(
                                                "Save Coordinates to Bed File",
                                                id="save-button",
                                                n_clicks=0,
                                                disabled=True,
                                                style={"margin": "20px"},
                                            ),
                                            dcc.Download(id="download"),
                                        ],
                                    ),
                                ],
                                style={
                                    "fontFamily": "Helvetica, Arial, sans-serif",
                                },
                            ),
                        ],
                        style={
                            "paddingTop": "110px",
                            "paddingLeft": "20px",
                            "textAlign": "left",
                        },  # Added paddingTop for the text and centered the text
                    ),
                ],
                style={
                    "display": "flex",
                    "justifyContent": "center",
                    "flex-wrap": "wrap",
                },
            ),
        ],
    )

    @app.callback(
        Output("text-input", "value"),
        Input("dotplot", "clickData"),
        Input("invisible-div", "children"),
        Input("matrix-dropdown", "value"),
    )
    def update_text(click_data, zoom_level, timothy):
        updated_info = ""
        for i in range(len(titles)):
            if timothy == titles[i]:
                updated_info = metadata[i]
        global clicked_values
        tmp = updated_info["max_window_size"] // (2**zoom_level)
        if click_data:
            x_val = (
                click_data["points"][0]["x"] if "x" in click_data["points"][0] else ""
            )
            y_val = (
                click_data["points"][0]["y"] if "y" in click_data["points"][0] else ""
            )
            z_val = (
                click_data["points"][0]["z"] if "z" in click_data["points"][0] else ""
            )

            # Create clicked_values list if it doesn't exist in the global scope
            if "clicked_values" not in globals():
                globals()["clicked_values"] = []

            # Append x and y values to the list and keep only the last 5 values
            # TODO: I shouldn't have to round, but I do :/
            clicked_values.insert(
                0,
                f'{updated_info["x_name"]}: {round(x_val)}-{round(x_val + tmp)}\n{updated_info["y_name"]}: {round(y_val)}-{round(y_val + tmp)}\nIdentity: {round(z_val,2)}\n~',
            )
            clicked_values = clicked_values[-100:]
            return "\n".join((clicked_values))
        else:
            return ""

    @app.callback(Output("save-button", "disabled"), Input("text-input", "value"))
    def update_button_state(text_value):
        return text_value == ""

    @app.callback(
        Output("text-input", "value", allow_duplicate=True),
        Output("save-bed", "n_clicks"),
        Input("save-bed", "n_clicks"),
        Input("dotplot", "figure"),
    )
    def save_bed(n_clicks, figure):
        global clicked_values
        if n_clicks > 0:
            window_size = figure["data"][0]["x"][1] - figure["data"][0]["x"][0]
            try:
                identity = round(
                    min([val for val in figure["data"][0]["z"][0] if val > 0])
                )
                x_axis_name = figure["layout"]["xaxis"]["title"]["text"]
                y_axis_name = figure["layout"]["yaxis"]["title"]["text"]
                if x_axis_name == y_axis_name:
                    selfy = True
                    title_hi = x_axis_name + ".bed"
                else:
                    selfy = False
                    title_hi = x_axis_name + "-" + y_axis_name + ".bed"
            except:
                identity = 86
                x_axis_name = "x"
                y_axis_name = "y"
                selfy = False
                title_hi = "x-y.bed"
            pls = np.array(figure["data"][0]["z"])
            tr = convertMatrixToBed(
                pls, window_size, identity, x_axis_name, y_axis_name, selfy
            )
            if not output_dir:
                bedfile_output = os.path.join("./", title_hi)
            else:
                bedfile_output = os.path.join(output_dir, title_hi)
                if (output_dir) and not os.path.exists(output_dir):
                    os.makedirs(output_dir)
            with open(bedfile_output, "w") as bedfile:
                for row in tr:
                    bedfile.write("\t".join(map(str, row)) + "\n")
            msg = f"Saved bed file to {bedfile_output}\n"
            return msg, 0  # Make sure to return a tuple of values for the outputs
        else:
            return (
                dash.no_update,
                dash.no_update,
            )  # Return dash.no_update for outputs not being updated

    @app.callback(
        Output("text-input", "value", allow_duplicate=True),
        Input("save-button", "n_clicks"),
        State("text-input", "value"),
    )
    def save_to_file(n_clicks, content):
        global clicked_values
        if n_clicks > 0 and content:
            content = content.replace("\n", ",")
            content = content.replace(": ", ", ")
            content = content.replace("~,", "\n")
            content = content.replace(" Identity,", "")
            content = content.replace("~", "")
            content = content.replace(",\n", "\n")
            msg = f"Saved coordinate history to coordinate_log.txt!"
            with open("coordinate_log.txt", "w") as file:
                file.write(str(content))
            return msg

    # Callback activates when panning, zooming, or changing color
    @app.callback(
        Output("dotplot", "figure", allow_duplicate=True),
        Output("loading-output-2", "children"),
        Output("window-div", "children"),
        Output("invisible-div", "children"),
        Input("dotplot", "relayoutData"),
        Input("dotplot", "figure"),
        Input("color-options", "value"),
        Input("threshold-slider", "value"),
        Input("gradient-toggle", "value"),
        Input("matrix-dropdown", "value"),
        Input("invisible-div", "children"),
        prevent_initial_call=True,
    )
    def update_dotplot(
        relayoutData,
        figure,
        color,
        threshold_range,
        button_states,
        timothy,
        current_zoom,
    ):
        updated_info = ""
        updated_title = ""
        title_size = 28
        image_axes = []
        # Matrix dropdown value = Timothy
        for i in range(len(titles)):
            if timothy == titles[i]:
                image_pyramid = matrices[i]
                image_axes = axes[i]
                updated_info = metadata[i]
                if updated_info["self"]:
                    updated_title = f"Self-Identity Plot: {updated_info['x_name']}"
                    if len(updated_info["x_name"]) > 22:
                        title_size = 20
                else:
                    updated_title = f"Comparative Plot: {updated_info['x_name']} vs. {updated_info['y_name']}"
                    if len(updated_info["x_name"]) + len(updated_info["y_name"]) > 22:
                        title_size = 16
        new_main_level = image_pyramid[0]
        fig = go.Figure(data=[heatmap])
        current_color = getInteractiveColor(getMatchingColors(color), "+")
        new_heatmap = heatmap
        new_heatmap.update(dict(colorscale=current_color, z=new_main_level))

        masked_data = np.where(
            (new_heatmap["z"] < threshold_range[0])
            | (new_heatmap["z"] > threshold_range[1]),
            0,
            new_heatmap["z"],
        )
        fig = go.Figure(data=[new_heatmap])
        fig.update_traces(z=masked_data)

        if "keep-original" in button_states:
            fig.update_traces(zmin=threshold_range[0])
            fig.update_traces(zmax=threshold_range[1])
        else:
            fig.update_traces(zmin=identity)
            fig.update_traces(zmax=100)

        fig.update_xaxes(
            showspikes=True,
            spikemode="across",
            ticks="outside",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
        )
        # TODO: nticks not working
        fig.update_xaxes(
            nticks=10, title_text=updated_info["x_name"], title_font=dict(size=18)
        )

        # Update y-axis properties directly within the heatmap trace
        fig.update_yaxes(
            showspikes=True,
            spikemode="across",
            ticks="outside",
            showline=True,
            linewidth=2,
            linecolor="black",
            mirror=True,
        )
        fig.update_yaxes(title_text=updated_info["y_name"], title_font=dict(size=18))

        fig.update_layout(
            hoverlabel=dict(bgcolor="white", font_size=16, font_family="Helvetica"),
            title=updated_title,  # Add your title here
            title_font=dict(
                size=title_size, family="Helvetica"
            ),  # Adjust the title font size if needed
            title_x=0.5,
            title_y=0.95,
        )

        fig.update_layout(yaxis_scaleanchor="x")
        if relayoutData is not None:
            # TODO: Pan mode should stay in pan mode

            if "xaxis.range[0]" in relayoutData:
                x_start_range = relayoutData.get("xaxis.range[0]")
                y_start_range = relayoutData.get("yaxis.range[0]")
                x_end_range = relayoutData.get("xaxis.range[1]")
                y_end_range = relayoutData.get("yaxis.range[1]")

                # Check that selected range is in bounds, snap to boundary otherwise
                # TODO: still broken here
                if x_start_range is not None:
                    if x_start_range < 0:
                        print("x axis out of bounds! Shifting to x=0")
                        relayoutData["xaxis.range[0]"] = 0
                        x_start_range = 0
                    if y_start_range < 0:
                        print("y axis out of bounds! Shifting to y=0")
                        relayoutData["yaxis.range[0]"] = 0
                        y_start_range = 0
                    if x_end_range > updated_info["x_size"]:
                        print(
                            f"x axis out of bounds! Shifting to x={updated_info['x_size']}"
                        )
                        relayoutData["xaxis.range[1]"] = updated_info["x_size"]
                        x_end_range = updated_info["x_size"]
                    if y_end_range > updated_info["y_size"]:
                        print(
                            f"y axis out of bounds! Shifting to y={updated_info['y_size']}"
                        )
                        relayoutData["yaxis.range[1]"] = updated_info["y_size"]
                        y_end_range = updated_info["y_size"]

                if x_start_range >= 0 and y_start_range >= 0:
                    x_begin = round(relayoutData["xaxis.range[0]"])
                    x_end = round(relayoutData["xaxis.range[1]"])

                    y_begin = round(relayoutData["yaxis.range[0]"])
                    y_end = round(relayoutData["yaxis.range[1]"])

                    amount = x_end - x_begin

                    # This function finds the correct level in the image pyramid
                    try:
                        zoom_factor = findValueInRange(amount, important)
                    except ValueError as err:
                        zoom_factor = 0
                    # If zoom_factor is less than current_zoom, base tmp factors on the older amount

                    if zoom_factor > len(image_pyramid) - 1:
                        zoom_factor = len(image_pyramid) - 1
                    if (zoom_factor != 0) or (zoom_factor != current_zoom):
                        using_matrix = image_pyramid[zoom_factor]
                        # Get the coordinates for the better matrix
                        tmp1 = find_closest_elements(
                            x_start_range, image_axes[2 * zoom_factor]
                        )[0]
                        tmp2 = find_closest_elements(
                            x_end_range, image_axes[2 * zoom_factor]
                        )[1]
                        tmp3 = find_closest_elements(
                            y_start_range, image_axes[(2 * zoom_factor) + 1]
                        )[0]
                        tmp4 = find_closest_elements(
                            y_end_range, image_axes[(2 * zoom_factor) + 1]
                        )[1]
                        x_start_updated = math.floor(tmp1)
                        x_end_updated = math.ceil(tmp2)
                        y_start_updated = math.floor(tmp3)
                        y_end_updated = math.ceil(tmp4)

                        x_equi, y_equi = makeDifferencesEqual(
                            x_start_updated,
                            x_end_updated,
                            y_start_updated,
                            y_end_updated,
                        )
                        # TODO: flip axes here
                        using_matrix = np.copy(
                            using_matrix[y_start_updated:y_equi, x_start_updated:x_equi]
                        )
                        masked_matrix = np.where(
                            (using_matrix < threshold_range[0])
                            | (using_matrix > threshold_range[1]),
                            0,
                            using_matrix,
                        )

                        x_ax = image_axes[2 * zoom_factor][x_start_updated:x_equi]
                        y_ax = image_axes[(2 * zoom_factor) + 1][y_start_updated:y_equi]

                        new_hover_template_text = "Identity: %{z:.1f}"

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
                            hoverinfo="z",  # Hover info to display
                            hovertemplate=new_hover_template_text,
                            name="",
                            x0=x_ax[0],
                            dx=updated_info["max_window_size"] // (2**zoom_factor),
                            xtype="scaled",
                            y0=y_ax[0],
                            dy=updated_info["max_window_size"] // (2**zoom_factor),
                            ytype="scaled",
                        )

                        # Create the figure
                        current_fig = go.Figure(data=[new_heatmap])
                        current_fig.update_xaxes(
                            showspikes=True,
                            spikemode="across",
                            ticks="outside",
                            showline=True,
                            linewidth=2,
                            linecolor="black",
                            mirror=True,
                        )
                        current_fig.update_xaxes(
                            nticks=10,
                            title_text=updated_info["x_name"],
                            title_font=dict(size=18),
                        )

                        # Update y-axis properties directly within the heatmap trace
                        current_fig.update_yaxes(
                            showspikes=True,
                            spikemode="across",
                            ticks="outside",
                            showline=True,
                            linewidth=2,
                            linecolor="black",
                            mirror=True,
                        )

                        current_fig.update_yaxes(
                            title_text=updated_info["y_name"], title_font=dict(size=18)
                        )
                        current_fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white", font_size=16, font_family="Helvetica"
                            ),
                            title=fig_title,  # Add your title here
                            title_font=dict(size=28, family="Helvetica"),
                            title_x=0.5,  # Center horizontally
                            title_y=0.95,  # Center vertically
                        )
                        current_fig.update_layout(yaxis_scaleanchor="x")
                        return (
                            current_fig,
                            "",
                            f"Current Window Size: {updated_info['max_window_size']//(2**zoom_factor)}",
                            zoom_factor,
                        )
                    else:
                        # Happens when panning or zooming in/out within threshold
                        return (
                            figure,
                            "",
                            f"Current Window Size: {updated_info['max_window_size']//(2**current_zoom)}",
                            current_zoom,
                        )

            elif "xaxis.autorange" in relayoutData:
                # Triggers when double-clicked or reset axes
                return (
                    fig,
                    "",
                    f"Current Window Size: {updated_info['max_window_size']}",
                    0,
                )

            # TODO: disable Autosizing
            elif "autosize" in relayoutData:
                # Autosize should be disabled: if not it gets handled here
                return (
                    figure,
                    "",
                    f"Current Window Size: {updated_info['max_window_size']}",
                    zoom_factor,
                )

            else:
                # Happens when user selects pan mode.
                # TODO: get correct window size
                return (
                    figure,
                    "",
                    f"Current Window Size: {updated_info['max_window_size']//(2**current_zoom)}",
                    current_zoom,
                )
        else:
            # This occurs at startup and when selecting color. Return the base figure with color
            # fig.update_layout(coloraxis=dict(colorscale=current_color))
            return (
                fig,
                "",
                f"Current Window Size: {updated_info['max_window_size']}",
                current_zoom,
            )

    # -------DO NOT DELETE! CUSTOM CSS FOR IDENTITY THRESHOLD SLIDER--------
    app.index_string = """
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
            #main_color {
                justify-content:center;
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
"""
    # TODO: Change debug to False for production
    app.run(debug=False, use_reloader=False, port=port_number)
