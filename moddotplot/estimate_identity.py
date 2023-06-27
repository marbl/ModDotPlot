#!/usr/bin/env python3
import math
import numpy as np
from moddotplot.const import (
    SEQUENTIAL_PALETTES,
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
)
from palettable import colorbrewer


def convert_set(kmer_list):
    return [set(k) for k in kmer_list]


# Update to just using the neighbors, not the whole set + neighbors!
def convert_set_neighbors(mod_list, alpha):
    tmp = []
    length = len(mod_list)
    for k in range(length):
        ll = set(mod_list[k])
        if k != 0:
            prev_length = round(len(mod_list[k - 1]) * alpha)
            ll.update(
                mod_list[k - 1][:prev_length]
            )  # Update the set with elements from the previous kmer_list
        if k != length - 1:
            high_length = round(len(mod_list[k - 1]) * alpha)
            ll.update(
                mod_list[k + 1][high_length:]
            )  # Update the set with elements from the next kmer_list
        tmp.append(ll)
    return tmp


def initial_identity_matrix(mod_set, resolution, k):
    n = len(mod_set)
    jaccard_matrix = np.empty((resolution, resolution))

    for w in range(n):
        for r in range(w, n):  # Iterate over the upper triangular matrix
            j_hat = binomial_distance(containment(mod_set[w], mod_set[r]), k)
            jaccard_matrix[w, r] = j_hat
            jaccard_matrix[r, w] = j_hat  # Symmetric assignment

    return jaccard_matrix


def updated_identity_matrix(
    mod_set_x, mod_set_y, mod_set_x_neighbors, mod_set_y_neighbors, resolution, k
):
    jaccard_matrix = np.empty((resolution, resolution))
    for w in range(len(mod_set_y)):
        for q in range(len(mod_set_x)):
            jaccard_matrix[w, q] = binomial_distance(
                containment_neighbors(
                    mod_set_x[q],
                    mod_set_y[w],
                    mod_set_x_neighbors[q],
                    mod_set_y_neighbors[w],
                ),
                k,
            )

    return jaccard_matrix


def find_elements_with_prefix(lst, prefix):
    matching_elements = []
    for element in lst:
        if element.startswith(prefix):
            matching_elements.append(element)
    return matching_elements


def get_interactive_color(palette_name, palette_orientation):
    palettes = colorbrewer.COLOR_MAPS
    tmp_color = []
    new_palette = palette_name.split("_")
    if palette_name in DIVERGING_PALETTES:
        tmp_color = palettes["Diverging"][new_palette[0]][new_palette[1]]["Colors"]
        if palette_orientation == "+":
            palette_orientation = "-"
        else:
            palette_orientation = "+"
    elif palette_name in SEQUENTIAL_PALETTES:
        tmp_color = palettes["Sequential"][new_palette[0]][new_palette[1]]["Colors"]
    elif palette_name in QUALITATIVE_PALETTES:
        tmp_color = palettes["Qualitative"][new_palette[0]][new_palette[1]]["Colors"]
    else:
        print("Unable to determine color palette. Selecting default \n")
        tmp_color = palettes["Diverging"]["Spectral"]["11"]["Colors"]
        palette_orientation = "-"
    if palette_orientation == "-":
        tmp_color = tmp_color[::-1]
    tmp_color = [[255, 255, 255]] + tmp_color
    total_values = len(tmp_color)
    formatted_values = [
        [i / (total_values - 1), f"rgb({r}, {g}, {b})"]
        for i, (r, g, b) in enumerate(tmp_color)
    ]
    return formatted_values


def get_matching_colors(color_name):
    available_colors = [
        element
        for sublist in [DIVERGING_PALETTES, QUALITATIVE_PALETTES, SEQUENTIAL_PALETTES]
        for element in sublist
    ]
    matching_elements = find_elements_with_prefix(available_colors, color_name)
    return matching_elements[-1]


def divide_into_chunks(lst, res):
    chunk_size = math.ceil(len(lst) / res)
    chunks = [lst[i : i + chunk_size] for i in range(0, len(lst), chunk_size)][:res]
    while len(chunks) > res:
        chunks.pop()
    while len(chunks) < res:
        chunks.append([])
    return chunks


def get_mods(kmer_list, s, res):
    assert int(s) > 0
    mod_list_prep = divide_into_chunks(kmer_list, res)
    mod_list = [
        [kmer for kmer in element if kmer % s == 0] for element in mod_list_prep
    ]
    return mod_list


def partition_evenly_spaced_modimizers(mod_list, seq_length, resolution):
    thousand_dict = {}
    for i in range(resolution):
        start_site = i * round(seq_length / resolution)
        end_site = (
            i * round(seq_length / resolution) + round(seq_length / resolution) - 1
        )
        name = f"{start_site}-{end_site}"
        thousand_dict[name] = mod_list[i]
    return thousand_dict


def create_coordinates(kmer_list, resolution):
    zoomed_limit = math.floor(len(kmer_list) / resolution)
    tmp = []
    for k in range(resolution):
        tmp.append(set())
        for j in range(zoomed_limit):
            tmp[k].add(kmer_list[(k * zoomed_limit) + j])

    return tmp


def create_coordinates2(kmer_list):
    tmp = []
    for k in kmer_list:
        tmp.append(set())
        for j in k:
            tmp[-1].add(j)

    return tmp


def binomial_distance(containment_value, kmer_value):
    return math.pow(containment_value, (1 / kmer_value))


def containment(set1, set2):
    intersection = set1.intersection(set2)
    try:
        if len(set1) > len(set2):
            return float(len(intersection) / len(set1))
        else:
            return float(len(intersection) / len(set2))
    except:
        return 0


def containment_neighbors(set1, set2, set3, set4):
    intersection = set3.intersection(set4)
    denominator = len(set1) if len(set1) > len(set2) else len(set2)
    try:
        result = float(len(intersection) / denominator)
        return 1 if result > 1 else result
    except ZeroDivisionError:
        return 0


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
    zoom_levels[sparsity_layers[0]] = axis_length
    for i in range(1, len(sparsity_layers)):
        zoom_levels[sparsity_layers[i]] = round(axis_length / pow(2, i))
    return zoom_levels


def jaccard(set1, set2):
    intersection = set1.intersection(set2)
    try:
        return float(len(intersection) / (len(set1) + len(set2) - len(intersection)))
    except:
        return 0


def poisson_distance(jaccard_value, kmer_value):
    return (-1 / kmer_value) * math.log((2 * jaccard_value / (1 + jaccard_value)))
