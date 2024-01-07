#!/usr/bin/env python3
import math
import numpy as np
from moddotplot.const import (
    SEQUENTIAL_PALETTES,
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
)
from palettable import colorbrewer
from typing import List, Set, Dict, Tuple

from moddotplot.parse_fasta import printProgressBar

def partition_overlaps(lst: List[int], win: int, delta: float, seq_len: int, k: int) -> List[List[int]]:

    kmer_list = []
    kmer_to_genomic_coordinate_offset = win - k + 1
    delta_offset = win * delta

    # Set the first window to contain win - k + 1 kmers.
    starting_end_index = int(round(kmer_to_genomic_coordinate_offset + delta_offset))
    kmer_list.append(lst[0:starting_end_index])
    counter = win - k + 1

    # Set normal windows 
    while counter <= (seq_len - win):
        start_index = counter + 1
        end_index = win + counter + 1
        delta_start_index = int(round(start_index - delta_offset))
        delta_end_index = int(round(end_index + delta_offset))
        if delta_end_index > seq_len:
            delta_end_index = seq_len
        try:
            kmer_list.append(lst[delta_start_index:delta_end_index])
        except Exception as e:
            print(e)
            kmer_list.append(lst[delta_start_index:seq_len])
        counter += win

    # Set the last window to get the remainder
    if (counter <= seq_len - 2):

        final_start_index = int(round(counter + 1 - delta_offset))
        kmer_list.append(lst[final_start_index:seq_len])

    # Test that last value was added on correctly

    assert kmer_list[-1][-1] == lst[-1]
    return kmer_list

def convert_to_modimizers(kmer_list: List[List[int]], sparsity: int) -> List[List[int]]:
    mod_tot = []
    for partition in kmer_list:
        mod_list = []
        for kmer in partition:
            if kmer % sparsity == 0:
                mod_list.append(kmer)
        #print(f"{len(mod_list)} -> {len(set(mod_list))}")
        mod_tot.append(set(mod_list))
    return mod_tot



def divide_into_chunks(lst: List[int], res: int) -> List[List[int]]:
    """
    Divide a list into approximately equal-sized chunks.

    Args:
        lst (List[int]): The input list to be divided.
        res (int): The desired number of result chunks.

    Returns:
        List[List[int]]: A list of lists, where each inner list contains elements from the input list.
    """
    chunk_size = len(lst) // res  # Calculate the target chunk size
    remainder = len(lst) % res  # Calculate the remainder

    chunks = []
    start = 0

    for i in range(res):
        end = (
            start + chunk_size + (1 if i < remainder else 0)
        )  # Adjust chunk size for remainder
        chunks.append(lst[start:end])
        start = end

    return chunks


def get_mods(kmer_list: List[int], s: int, res: int) -> List[List[int]]:
    """
    Filters a list of k-mers based on divisibility by s.

    Args:
        kmer_list (List[int]): A list of integers representing k-mers.
        s (int): A positive integer used as a divisor.
        res (int): The desired number of result chunks.

    Returns:
        List[List[int]]: A list of lists, where each inner list contains k-mers from kmer_list
                         that are divisible by s.
    """
    if s < 1:
        s = 1

    mod_list_prep = divide_into_chunks(kmer_list, res)

    mod_list = [
        [kmer for kmer in element if kmer % s == 0] for element in mod_list_prep
    ]

    return mod_list


def convert_set(kmer_list: List[List[int]]) -> List[Set[int]]:
    """
    Convert a list of lists of integers into a list of sets.

    Args:
        kmer_list (List[List[int]]): A list of lists, where each inner list contains integers.

    Returns:
        List[Set[int]]: A list of sets, where each set contains integers from the input lists.
    """
    return [set(k) for k in kmer_list]


def convert_set_neighbors(mod_list: List[List[int]], delta: float) -> List[Set[int]]:
    """
    Convert a list of lists of integers into a list of sets with neighboring elements.

    Args:
        mod_list (List[List[int]]): A list of lists, where each inner list contains integers.
        delta (float): A scaling factor to determine the neighborhood size.

    Returns:
        List[Set[int]]: A list of sets, where each set contains integers from neighboring elements in mod_list.
    """
    kmer_sets = []
    length = len(mod_list)

    for k in range(length):
        ll = set(mod_list[k])

        if k != 0:
            prev_length = round(len(mod_list[k - 1]) * delta)
            print("\nprev")
            print(prev_length)
            ll.update(mod_list[k - 1][-prev_length:])

        if k != length - 1:
            high_length = round(len(mod_list[k + 1]) * delta)
            print("\nhigh")
            print(high_length)
            ll.update(mod_list[k + 1][:high_length])
        
        kmer_sets.append(ll)

    return kmer_sets


def binomial_distance(containment_value: float, kmer_value: int) -> float:
    """
    Calculate the binomial distance based on containment and kmer values.

    Args:
        containment_value (float): The containment value.
        kmer_value (int): The k-mer value.

    Returns:
        float: The binomial distance.
    """
    return math.pow(containment_value, 1.0 / kmer_value)


def containment_neighbors(
    set1: Set[int],
    set2: Set[int],
    set3: Set[int],
    set4: Set[int],
    identity: int,
    k: int,
) -> float:
    """
    Calculate the containment neighbors based on four sets and an identity threshold.

    Args:
        set1 (Set[int]): The first set.
        set2 (Set[int]): The second set.
        set3 (Set[int]): The third set.
        set4 (Set[int]): The fourth set.
        identity (int): The identity threshold.
        k (int): Kmer value.

    Returns:
        float: The containment neighbors value.
    """
    len_a = len(set1)
    len_b = len(set2)

    intersection_a_b_prime = len(set1 & set4)
    if len_a != 0:
        containment_a_b_prime = intersection_a_b_prime / len_a
    else:
        # If len_a is zero, handle it by setting containment_a_b_prime to a default value
        containment_a_b_prime = 0

    if binomial_distance(containment_a_b_prime, k) < identity / 100:
        return 0.0

    else:
        intersection_a_prime_b = len(set2 & set3)
        if len_b != 0:
            containment_a_prime_b = intersection_a_prime_b / len_b
        else:
            # If len_a is zero, handle it by setting containment_a_b_prime to a default value
            containment_a_prime_b = 0

        return max(containment_a_b_prime, containment_a_prime_b)


def self_containment_matrix(
    mod_set: List[set], mod_set_neighbors: List[set], k: int, identity: int
) -> np.ndarray:
    """
    Create a self-containment matrix based on containment similarity calculations.

    Args:
        mod_set (List[set]): A list of sets representing elements.
        mod_set_neighbors (List[set]): A list of sets representing neighbors for each element.
        k (int): A parameter for containment similarity calculation.

    Returns:
        np.ndarray: A NumPy array representing the self-containment matrix.
    """
    n = len(mod_set)
    progress_thresholds = round(n / 77)

    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    containment_matrix = np.empty((n, n))

    for w in range(n):
        if w % progress_thresholds == 0:
            printProgressBar(w, n, prefix="Progress:", suffix="Complete", length=40)
        containment_matrix[w, w] = 1.0

        for r in range(w + 1, n):
            c_hat = binomial_distance(
                containment_neighbors(
                    mod_set[w],
                    mod_set[r],
                    mod_set_neighbors[w],
                    mod_set_neighbors[r],
                    identity,
                    k,
                ),
                k,
            )
            containment_matrix[r, w] = c_hat
            containment_matrix[w, r] = c_hat

    printProgressBar(
        n, n, prefix="Progress:", suffix="Completed", length=40
    )  # show completed progress bar
    print("\n")
    return containment_matrix


def pairwise_containment_matrix(
    mod_set_x: List[int],
    mod_set_y: List[int],
    mod_set_x_neighbors: List[List[int]],
    mod_set_y_neighbors: List[List[int]],
    identity: int,
    k: int,
    supress_progress: bool,
) -> np.ndarray:
    """
    Calculate an updated identity matrix using specified parameters.

    Args:
        mod_set_x (List[int]): List of values for the x-axis.
        mod_set_y (List[int]): List of values for the y-axis.
        mod_set_x_neighbors (List[List[int]]): List of lists representing neighbors of mod_set_x values.
        mod_set_y_neighbors (List[List[int]]): List of lists representing neighbors of mod_set_y values.
        identity (int): Resolution parameter.
        k (int): Value for the k parameter in the binomial_distance function.
        supress_progress (bool): if true supresses the progress bar

    Returns:
        np.ndarray: An identity matrix containing containment values.
    """
    #X is the larger, y is the smaller
    print("sanity check")
    n = len(mod_set_x)
    progress_thresholds = round(n / 77)

    if not supress_progress:
        printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    containment_matrix = np.zeros((n, n), dtype=float)

    for w in range(len(mod_set_y)):
        if not supress_progress:
            if w % progress_thresholds == 0:
                printProgressBar(w, n, prefix="Progress:", suffix="Complete", length=40)
        for q in range(n):
            containment_matrix[w, q] = binomial_distance(
                containment_neighbors(
                    mod_set_x[q],
                    mod_set_y[w],
                    mod_set_x_neighbors[q],
                    mod_set_y_neighbors[w],
                    identity,
                    k,
                ),
                k,
            )
    
    if not supress_progress:
        printProgressBar(
            n, n, prefix="Progress:", suffix="Completed", length=40
        )  # show completed progress bar
        print("\n")
    return containment_matrix


# Function used to find matching color palette to those available in const.py
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

def partition_evenly_spaced_modimizers(mod_list, seq_length, resolution):
    coords_dict = {}
    for i in range(resolution):
        start_site = i * round(seq_length / resolution) + 1
        end_site = (
            i * round(seq_length / resolution) + round(seq_length / resolution) - 1
        ) + 1
        name = f"{start_site}-{end_site}"
        coords_dict[name] = mod_list[i]
    return coords_dict

def partition_pairwise_modimizers_different_size(mod_list_large, mod_list_small, ratio, window_size):
    short_seq_dict = {}
    large_seq_dict = {}
    assert ratio == len(mod_list_small)
    for i in range(ratio):
        start_site = (i * window_size) + 1
        end_site = start_site + window_size - 1
        name = f"{start_site}-{end_site}"
        short_seq_dict[name] = mod_list_small[i]
        large_seq_dict[name] = mod_list_large[i]
    for j in range(len(mod_list_large) - ratio):
        start_site = (ratio * window_size) + (j * window_size) + 1
        end_site = start_site + window_size - 1
        name = f"{start_site}-{end_site}"
        large_seq_dict[name] = mod_list_large[j]
    return short_seq_dict, large_seq_dict



def containment(set1, set2):
    intersection = set1.intersection(set2)
    try:
        if len(set1) > 0 and len(set2) > 0:
            if len(set1) > len(set2):
                return float(len(intersection) / len(set1))
            else:
                return float(len(intersection) / len(set2))
        else:
            return 0.0
    except ZeroDivisionError:
        return 0.0


## ALL USEFUL


def verify_modimizers(sparsity, l):
    # Get the next highest power of 2, if not provided
    updated_sparsity = next_power_of_two(sparsity)

    sparsity_layers = [updated_sparsity]
    while l > 0:
        if sparsity_layers[-1] == 1:
            return sparsity_layers
        elif sparsity_layers[-1] % 2 == 1:
            sparsity_layers[-1] = int(sparsity_layers[-1] + 1)
        sparsity_layers.append(int(sparsity_layers[-1] / 2))
        l -= 1

    return sparsity_layers


def generate_dict_from_list(lst: List[int]) -> Dict[Tuple[int, int], int]:
    result = {}
    for i in range(len(lst) - 1):
        result[(lst[i], lst[i + 1])] = i
    return result


def find_value_in_range(integer: int, range_dict: dict) -> int:
    if integer > max(key[0] for key in range_dict.keys()):
        return 0
    highest_value = max(range_dict.values()) + 1
    for key, value in range_dict.items():
        if key[0] >= integer >= key[1]:
            return value
    return highest_value


## ALL USEFUL


def set_zoom_levels(axis_length, sparsity_layers):
    zoom_levels = {}
    zoom_levels[sparsity_layers[0]] = axis_length
    for i in range(1, len(sparsity_layers)):
        zoom_levels[sparsity_layers[i]] = round(axis_length / pow(2, i))
    return zoom_levels


def set_zoom_levels_list(axis_length, sparsity_layers):
    zoom_levels = []
    zoom_levels.append(axis_length)
    for i in range(1, len(sparsity_layers)):
        zoom_levels.append(round(axis_length / pow(2, i)))
    return zoom_levels


def find_closest_higher_index(descending_list, target_value):
    left = 0
    right = len(descending_list) - 1
    closest_index = None

    while left <= right:
        mid = (left + right) // 2

        if descending_list[mid] > target_value:
            closest_index = mid
            left = mid + 1
        else:
            right = mid - 1

    return closest_index


def make_differences_equal(x, x_prime, y, y_prime):
    difference_x = abs(x_prime - x)
    difference_y = abs(y_prime - y)

    if difference_x != difference_y:
        if difference_x < difference_y:
            x_prime += difference_y - difference_x
        else:
            y_prime += difference_x - difference_y

    return x_prime, y_prime


def next_power_of_two(n):
    if n <= 0:
        return 1
    n -= 1
    n |= n >> 1
    n |= n >> 2
    n |= n >> 4
    n |= n >> 8
    n |= n >> 16
    return n + 1
