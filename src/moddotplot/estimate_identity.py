#!/usr/bin/env python3
import math
from collections import OrderedDict
from dataclasses import dataclass
import numpy as np
from moddotplot.const import (
    SEQUENTIAL_PALETTES,
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
)
from palettable import colorbrewer
from typing import Collection, Hashable, List, Set, Dict, Tuple
import pandas as pd
import cooler
from scipy.sparse import csr_matrix

from moddotplot.parse_fasta import printProgressBar


@dataclass(frozen=True)
class PreparedModimizerSketches:
    """Core and expanded window sketches prepared for a matrix calculation.

    The compact sorted NumPy arrays are deliberately retained by reference. A
    prepared value can therefore be shared by self and pairwise calculations
    without copying sketch data. Public conversion helpers continue to return
    sets for backward compatibility.
    """

    core: List[Collection[int]]
    neighbors: List[Collection[int]]


def prepare_modimizer_sketches(
    sequence_length,
    sequence,
    window_size,
    sparsity,
    delta,
    k,
    ambiguous,
    expectation,
):
    """Partition one sequence and build the sketches used by matrix routines."""

    core_partitions = partitionOverlaps(sequence, window_size, 0, sequence_length, k)
    if delta > 0:
        neighbor_partitions = partitionOverlaps(
            sequence, window_size, delta, sequence_length, k
        )
    else:
        neighbor_partitions = core_partitions

    # Prepared values stay as compact sorted ndarrays.  The public conversion
    # helpers still return sets for API compatibility, but keeping millions of
    # selected hashes as Python ints in Python hash tables costs roughly ten
    # times more memory on chromosome-sized inputs.
    core = _convert_to_modimizer_arrays(
        core_partitions, sparsity, ambiguous, k, expectation
    )
    if neighbor_partitions is core_partitions:
        neighbors = core
    else:
        neighbors = _convert_to_modimizer_arrays(
            neighbor_partitions, sparsity, ambiguous, k, expectation
        )
    return PreparedModimizerSketches(core=core, neighbors=neighbors)


def create_self_matrix_from_sketches(prepared, k, identity, ambiguous):
    """Build a self matrix from an already prepared sequence sketch."""

    return selfContainmentMatrix(
        prepared.core, prepared.neighbors, k, identity, ambiguous
    )


def create_pairwise_matrix_from_sketches(
    prepared_x, prepared_y, identity, k, supress_progress=False
):
    """Build a pairwise matrix from two already prepared sequence sketches."""

    return pairwiseContainmentMatrix(
        prepared_x.core,
        prepared_y.core,
        prepared_x.neighbors,
        prepared_y.neighbors,
        identity,
        k,
        supress_progress,
    )


class ModimizerSketchCache:
    """Small LRU cache for prepared sequence sketches.

    ``source_key`` identifies the underlying sequence slice.  Calculation
    parameters are incorporated automatically, preventing reuse when a
    pairwise plot uses a different window size from the corresponding self
    plot. The default capacity matches a two-sequence comparison and bounds
    retained memory for larger grids.
    """

    def __init__(self, max_entries=2):
        if max_entries < 1:
            raise ValueError("max_entries must be at least one")
        self.max_entries = max_entries
        self._cache = OrderedDict()

    def get_or_prepare(
        self,
        source_key: Hashable,
        sequence_length,
        sequence,
        window_size,
        sparsity,
        delta,
        k,
        ambiguous,
        expectation,
    ):
        cache_key = (
            source_key,
            sequence_length,
            window_size,
            sparsity,
            delta,
            k,
            ambiguous,
            expectation,
        )
        try:
            prepared = self._cache.pop(cache_key)
        except KeyError:
            # Evict before constructing the replacement so a miss cannot
            # briefly exceed the configured memory bound by one full sketch.
            if len(self._cache) >= self.max_entries:
                self._cache.popitem(last=False)
            prepared = prepare_modimizer_sketches(
                sequence_length,
                sequence,
                window_size,
                sparsity,
                delta,
                k,
                ambiguous,
                expectation,
            )
        self._cache[cache_key] = prepared
        return prepared

    def clear(self):
        """Release references to all retained sketches."""

        self._cache.clear()


def createSelfMatrix(
    sequence_length,
    sequence,
    window_size,
    sparsity,
    delta,
    k,
    identity,
    ambiguous,
    sketch_size,
):
    prepared = prepare_modimizer_sketches(
        sequence_length,
        sequence,
        window_size,
        sparsity,
        delta,
        k,
        ambiguous,
        sketch_size,
    )
    return create_self_matrix_from_sketches(prepared, k, identity, ambiguous)


def createPairwiseMatrix(
    larger_length,
    smaller_length,
    larger_seq,
    smaller_seq,
    window_size,
    sparsity,
    delta,
    k,
    identity,
    ambiguous,
    expectation,
):
    prepared_large = prepare_modimizer_sketches(
        larger_length,
        larger_seq,
        window_size,
        sparsity,
        delta,
        k,
        ambiguous,
        expectation,
    )
    prepared_small = prepare_modimizer_sketches(
        smaller_length,
        smaller_seq,
        window_size,
        sparsity,
        delta,
        k,
        ambiguous,
        expectation,
    )
    return create_pairwise_matrix_from_sketches(
        prepared_large, prepared_small, identity, k
    )


def partitionOverlaps(
    lst: List[int], win: int, delta: float, seq_len: int, k: int
) -> List[List[int]]:
    if win <= 0:
        raise ValueError("window size must be greater than zero")
    if k <= 0:
        raise ValueError("k-mer size must be greater than zero")

    kmer_count = min(len(lst), max(seq_len, 0))
    if kmer_count == 0:
        return []

    # A sequence with n bases has n - k + 1 k-mers.  Reconstruct the genomic
    # length so that every partition starts on a multiple of ``win``.  The old
    # counter started the second partition at ``win - k + 2`` and then advanced
    # by ``win``, making all but the first window begin k - 2 bases too early.
    sequence_length = kmer_count + k - 1
    delta_offset = win * delta
    kmer_list = []

    # A trailing genomic fragment shorter than k has no k-mer and therefore no
    # matrix cell, so iterate over valid k-mer starts rather than base length.
    for window_start in range(0, kmer_count, win):
        window_end = min(window_start + win, sequence_length)
        expanded_start = max(0, int(round(window_start - delta_offset)))
        expanded_end = min(sequence_length, int(round(window_end + delta_offset)))

        # K-mers are indexed by their genomic start.  Subtracting k - 1 from
        # the right boundary excludes k-mers that cross out of the interval.
        start_index = min(expanded_start, kmer_count)
        end_index = min(max(expanded_end - k + 1, start_index), kmer_count)
        kmer_list.append(lst[start_index:end_index])

    return kmer_list


def _valid_hashes(partition):
    """Return the unmasked hashes in *partition* as a flat NumPy array.

    FASTA input uses a ``uint64`` ``MaskedArray``.  Keeping that representation
    here is important: iterating over a masked array produces a Python object
    for every genomic k-mer and was the dominant cost of sketch construction.
    The object-array fallback retains compatibility with legacy callers that
    pass lists containing ``None`` or masked scalar values.
    """

    if isinstance(partition, np.ma.MaskedArray):
        values = np.asarray(partition.data).reshape(-1)
        raw_mask = partition.mask
        # ``getmaskarray`` materializes a full all-False array for ``nomask``.
        # Avoid allocating and scanning that temporary for every 100 kb window.
        if raw_mask is np.ma.nomask or (np.ndim(raw_mask) == 0 and not bool(raw_mask)):
            return values
        mask = np.asarray(raw_mask, dtype=bool).reshape(-1)
        if np.any(mask):
            values = values[~mask]
        return values

    values = np.asarray(partition)
    # NumPy may infer ``float64`` for a Python list containing uint64-range
    # integers on some supported versions, irreversibly rounding hash values.
    # Send non-integral legacy sequences through the exact Python-int fallback
    # below; native numeric arrays can still be converted in bulk.
    legacy_requires_exact_conversion = (
        not isinstance(partition, np.ndarray) and values.dtype.kind not in "biu"
    )
    if values.dtype.kind != "O" and not legacy_requires_exact_conversion:
        values = values.reshape(-1)
        # ``populateModimizers`` historically called int() on every value.
        # Hash arrays are already integral, but retain that behavior for older
        # callers that provide a floating-point sequence.
        if values.dtype.kind not in "biu":
            values = values.astype(np.int64)
        return values

    # Mixed Python sequences cannot be converted to a numeric array until
    # ``None`` and masked sentinels have been removed.  This path is for API
    # compatibility; normal FASTA processing always takes the vectorized path
    # above.
    object_values = np.asarray(partition, dtype=object).reshape(-1)
    valid_values = [
        int(kmer)
        for kmer in object_values
        if kmer is not None and not np.ma.is_masked(kmer)
    ]
    if not valid_values:
        return np.empty(0, dtype=np.uint64)

    try:
        return np.asarray(valid_values, dtype=np.uint64)
    except (OverflowError, ValueError):
        # Preserve arbitrary-size and negative Python integers for legacy
        # callers. NumPy still performs the modulo and uniqueness operations
        # in bulk, albeit with object arithmetic.
        return np.asarray(valid_values, dtype=object)


def _divisible_hashes(values, sparsity):
    """Select hashes divisible by an integer sparsity without Python loops."""

    if values.size == 0:
        return values

    # All production sparsities are powers of two.  A bit mask avoids creating
    # the full-size temporary remainder array and is valid for signed and
    # unsigned integer hashes alike.
    if values.dtype.kind in "biu" and sparsity & (sparsity - 1) == 0:
        return values[np.bitwise_and(values, sparsity - 1) == 0]
    return values[np.remainder(values, sparsity) == 0]


def _populate_modimizer_array(partition, sparsity, ambiguous, expectation, k):
    """Build one adaptive modimizer sketch using vectorized NumPy operations.

    ``ambiguous`` and ``k`` remain in the public signature for compatibility;
    ambiguity is represented by the mask on ``partition``.  If a sketch is
    smaller than half its expectation, sparsity is repeatedly halved just as
    before.  The sparsity now stays integer throughout that fallback rather
    than becoming a float after the first recursive call.
    """

    del ambiguous, k

    current_sparsity = int(sparsity)
    if current_sparsity < 1:
        raise ValueError("sparsity must be a positive integer")

    values = _valid_hashes(partition)
    minimum_size = round(expectation / 2)

    while True:
        selected = _divisible_hashes(values, current_sparsity)
        unique = np.unique(selected)
        if unique.size >= minimum_size or current_sparsity == 1:
            return unique
        current_sparsity = max(1, current_sparsity // 2)


def populateModimizers(partition, sparsity, ambiguous, expectation, k):
    """Build one adaptive sketch and return its historical ``set`` type."""

    unique = _populate_modimizer_array(partition, sparsity, ambiguous, expectation, k)
    # ``tolist`` converts NumPy integer scalars to Python ints, preserving the
    # public return type while internal prepared sketches retain compact arrays.
    return set(unique.tolist())


def _convert_to_modimizer_arrays(
    kmer_list, sparsity: int, ambiguous: bool, k: int, expectation: int
):
    return [
        _populate_modimizer_array(partition, sparsity, ambiguous, expectation, k)
        for partition in kmer_list
    ]


def convertToModimizers(
    kmer_list: List[List[int]], sparsity: int, ambiguous: bool, k: int, expectation: int
) -> List[Set[int]]:
    return [
        populateModimizers(partition, sparsity, ambiguous, expectation, k)
        for partition in kmer_list
    ]


def convertMatrixToBed(
    matrix,
    window_size,
    id_threshold,
    x_name,
    y_name,
    self_identity,
    x_offset,
    y_offset,
    x_end=None,
    y_end=None,
):
    bed = [
        (
            "#query_name",
            "query_start",
            "query_end",
            "reference_name",
            "reference_start",
            "reference_end",
            "perID_by_events",
        )
    ]

    rows, cols = matrix.shape
    for x in range(rows):
        for y in range(cols):
            value = matrix[x, y]
            if (not self_identity) or (self_identity and x <= y):
                if value >= id_threshold / 100:
                    start_x = x * window_size + x_offset
                    end_x = start_x + window_size - 1
                    start_y = y * window_size + y_offset
                    end_y = start_y + window_size - 1

                    # The final matrix window can be shorter than ``window_size``.
                    # Keep exported coordinates inside the exact sequence or
                    # requested region instead of allowing the last tile to
                    # overhang it.
                    if x_end is not None:
                        end_x = min(end_x, x_end)
                    if y_end is not None:
                        end_y = min(end_y, y_end)

                    bed.append(
                        (
                            x_name,
                            int(start_x),
                            int(end_x),
                            y_name,
                            int(start_y),
                            int(end_y),
                            float(value),
                        )
                    )
    return bed


def convertMatrixToCool(
    matrix,
    window_size,
    id_threshold,
    x_name,
    y_name,
    self_identity,
    x_offset,
    y_offset,
    chromsizes,
    output_cool,
):
    """
    Convert a matrix into a .cool file.

    Args:
        matrix (ndarray): 2D numpy array.
        window_size (int): Bin/window size in bp.
        id_threshold (float): Percent identity threshold (0-100).
        x_name (str): Chromosome name for rows.
        y_name (str): Chromosome name for columns.
        self_identity (bool): Whether to include self and upper triangle only.
        x_offset (int): Genomic offset for x axis (start position).
        y_offset (int): Genomic offset for y axis (start position).
        chromsizes (dict): Dict of chromosome lengths, e.g. {"chr1": 248956422}.
        output_cool (str): Path to save cooler file.
    """
    rows, cols = matrix.shape

    # ---- build bin table ----
    bins = []
    for i in range(rows):
        bins.append(
            (x_name, i * window_size + x_offset, (i + 1) * window_size + x_offset)
        )
    for j in range(cols):
        bins.append(
            (y_name, j * window_size + y_offset, (j + 1) * window_size + y_offset)
        )

    bins = pd.DataFrame(bins, columns=["chrom", "start", "end"])

    # ---- build pixel table ----
    pixels = []
    for x in range(rows):
        for y in range(cols):
            value = matrix[x, y]
            if (not self_identity) or (self_identity and x <= y):
                if value >= id_threshold / 100:
                    bin1_id = x
                    bin2_id = rows + y  # offset y bins after x bins
                    pixels.append((bin1_id, bin2_id, float(value)))

    pixels = pd.DataFrame(pixels, columns=["bin1_id", "bin2_id", "count"])

    # ---- write cooler ----
    cooler.create_cooler(output_cool, bins=bins, pixels=pixels, ordered=True)
    print(bins)
    print(pixels)

    return output_cool


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
    Calculate symmetric containment using the opposite expanded window.

    ``set1`` and ``set2`` are sketches of the two core windows, while ``set3``
    and ``set4`` are the corresponding sketches expanded by ``delta``.  A core
    is compared with the *other* window's expanded sketch in both directions.
    This is what lets a repeat that straddles a partition boundary match a core
    window instead of being missed solely because the partitions are offset.

    The maximum directional containment is thresholded after both directions
    have been evaluated.  Applying the cutoff to only the first direction made
    the result depend on argument order in the original implementation.

    Args:
        set1 (Set[int]): The first set.
        set2 (Set[int]): The second set.
        set3 (Set[int]): Expanded sketch corresponding to ``set1``.
        set4 (Set[int]): Expanded sketch corresponding to ``set2``.
        identity (int): The identity threshold.
        k (int): Kmer value.

    Returns:
        float: The containment neighbors value.
    """
    containment_a_b_expanded = len(set1 & set4) / len(set1) if set1 else 0.0
    containment_b_a_expanded = len(set2 & set3) / len(set2) if set2 else 0.0
    symmetric_containment = max(containment_a_b_expanded, containment_b_a_expanded)

    if binomial_distance(symmetric_containment, k) < identity / 100:
        return 0.0
    return symmetric_containment


def _sketch_intersection_counts(sketches_a, sketches_b):
    """Return every pairwise sketch intersection count as a dense array.

    The previous matrix implementation performed two Python set intersections
    for every output cell.  At the default resolution that means roughly two
    million intersections per matrix, each scanning about 1,600 hashes.  Here
    hashes are coordinate-compressed once and the sketches are represented as
    sparse incidence matrices.  Sparse matrix multiplication then calculates
    the exact same intersection counts in compiled code.

    Only the small ``len(sketches_a) x len(sketches_b)`` result is dense.  The
    incidence matrices retain one entry per selected hash, so chromosome size
    does not create a dense hash universe.  Production ntHash values take the
    fast unsigned-64-bit path; the mapping fallback preserves exact behavior
    for legacy callers using negative, oversized, or other hashable values.
    """

    rows = len(sketches_a)
    cols = len(sketches_b)
    lengths_a = np.fromiter(
        (len(sketch) for sketch in sketches_a), dtype=np.int64, count=rows
    )
    lengths_b = np.fromiter(
        (len(sketch) for sketch in sketches_b), dtype=np.int64, count=cols
    )
    total_a = int(lengths_a.sum())
    total_b = int(lengths_b.sum())
    total = total_a + total_b

    count_dtype = np.int32
    if total == 0:
        return np.zeros((rows, cols), dtype=count_dtype)

    collections = (sketches_a, sketches_b)
    compact_uint64_arrays = all(
        isinstance(sketch, np.ndarray)
        and sketch.ndim == 1
        and sketch.dtype == np.dtype(np.uint64)
        for sketches in collections
        for sketch in sketches
    )

    if compact_uint64_arrays:
        # This is the normal prepared-sketch path. Concatenating the arrays in
        # C avoids boxing millions of hashes back into Python integers.
        hashes = np.concatenate(
            tuple(sketch for sketches in collections for sketch in sketches)
        )
        unique_hashes, inverse = np.unique(hashes, return_inverse=True)
        hash_count = len(unique_hashes)
        del hashes, unique_hashes
    else:
        uint64_max = np.iinfo(np.uint64).max
        uint64_compatible = all(
            isinstance(value, (int, np.integer)) and 0 <= int(value) <= uint64_max
            for sketches in collections
            for sketch in sketches
            for value in sketch
        )

    if not compact_uint64_arrays and uint64_compatible:
        hashes = np.fromiter(
            (
                int(value)
                for sketches in collections
                for sketch in sketches
                for value in sketch
            ),
            dtype=np.uint64,
            count=total,
        )
        unique_hashes, inverse = np.unique(hashes, return_inverse=True)
        hash_count = len(unique_hashes)
        del hashes, unique_hashes
    elif not compact_uint64_arrays:
        # This compatibility path is not used by FASTA processing, but keeps
        # the public matrix functions exact for arbitrary hashable set values.
        hash_columns = {}
        inverse = np.empty(total, dtype=np.int64)
        position = 0
        for sketches in collections:
            for sketch in sketches:
                for value in sketch:
                    try:
                        column = hash_columns[value]
                    except KeyError:
                        column = len(hash_columns)
                        hash_columns[value] = column
                    inverse[position] = column
                    position += 1
        hash_count = len(hash_columns)

    # SciPy uses 32-bit sparse indices whenever dimensions and nonzero counts
    # fit.  Keeping that representation saves tens of megabytes at r=1000.
    index_dtype = (
        np.int32
        if hash_count <= np.iinfo(np.int32).max and total <= np.iinfo(np.int32).max
        else np.int64
    )
    indices_a = inverse[:total_a].astype(index_dtype, copy=False)
    indices_b = inverse[total_a:].astype(index_dtype, copy=False)
    del inverse

    indptr_a = np.empty(rows + 1, dtype=index_dtype)
    indptr_b = np.empty(cols + 1, dtype=index_dtype)
    indptr_a[0] = 0
    indptr_b[0] = 0
    np.cumsum(lengths_a, out=indptr_a[1:])
    np.cumsum(lengths_b, out=indptr_b[1:])

    incidence_a = csr_matrix(
        (
            np.ones(total_a, dtype=count_dtype),
            indices_a,
            indptr_a,
        ),
        shape=(rows, hash_count),
    )
    incidence_b = csr_matrix(
        (
            np.ones(total_b, dtype=count_dtype),
            indices_b,
            indptr_b,
        ),
        shape=(cols, hash_count),
    )
    return (incidence_a @ incidence_b.T).toarray()


def _identity_matrix_from_containment(containment_matrix, identity, k):
    """Apply ModDotPlot's k-mer identity transform and cutoff in place."""

    identities = np.power(containment_matrix, 1.0 / k)
    identities[identities < identity / 100] = 0.0
    identities *= 100.0
    return identities


def selfContainmentMatrix(
    mod_set: List[set],
    mod_set_neighbors: List[set],
    k: int,
    identity: int,
    ambiguous: bool,
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
    if len(mod_set_neighbors) != n:
        raise IndexError("core and expanded self sketches must have equal lengths")
    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    intersection_counts = _sketch_intersection_counts(mod_set, mod_set_neighbors)
    core_sizes = np.fromiter((len(sketch) for sketch in mod_set), dtype=float, count=n)
    directional_containment = np.zeros((n, n), dtype=float)
    np.divide(
        intersection_counts,
        core_sizes[:, np.newaxis],
        out=directional_containment,
        where=core_sizes[:, np.newaxis] != 0,
    )
    symmetric_containment = np.maximum(
        directional_containment, directional_containment.T
    )
    containment_matrix = _identity_matrix_from_containment(
        symmetric_containment, identity, k
    )

    diagonal = np.full(n, 100.0)
    if not ambiguous:
        diagonal[core_sizes == 0] = 0.0
    np.fill_diagonal(containment_matrix, diagonal)

    printProgressBar(
        n, n, prefix="Progress:", suffix="Completed", length=40
    )  # show completed progress bar
    print("\n")
    return containment_matrix


def pairwiseContainmentMatrix(
    mod_set_x: List[Set[int]],
    mod_set_y: List[Set[int]],
    mod_set_x_neighbors: List[Set[int]],
    mod_set_y_neighbors: List[Set[int]],
    identity: int,
    k: int,
    supress_progress: bool,
) -> np.ndarray:
    """
    Calculate an updated identity matrix using specified parameters.

    Args:
        mod_set_x (List[Set[int]]): Modimizer sets for columns on the x-axis.
        mod_set_y (List[Set[int]]): Modimizer sets for rows on the y-axis.
        mod_set_x_neighbors (List[Set[int]]): Neighbor sets for x-axis windows.
        mod_set_y_neighbors (List[Set[int]]): Neighbor sets for y-axis windows.
        identity (int): Resolution parameter.
        k (int): Value for the k parameter in the binomial_distance function.
        supress_progress (bool): if true supresses the progress bar

    Returns:
        np.ndarray: A ``(len(mod_set_y), len(mod_set_x))`` identity matrix.
    """
    rows = len(mod_set_y)
    cols = len(mod_set_x)
    if len(mod_set_x_neighbors) != cols or len(mod_set_y_neighbors) != rows:
        raise IndexError("core and expanded pairwise sketches must have equal lengths")
    if not supress_progress:
        printProgressBar(0, rows, prefix="Progress:", suffix="Complete", length=40)
    x_core_sizes = np.fromiter(
        (len(sketch) for sketch in mod_set_x), dtype=float, count=cols
    )
    y_core_sizes = np.fromiter(
        (len(sketch) for sketch in mod_set_y), dtype=float, count=rows
    )

    # core X against expanded Y, transposed into the public (Y, X) layout.
    x_to_y_counts = _sketch_intersection_counts(mod_set_x, mod_set_y_neighbors).T
    x_to_y = np.zeros((rows, cols), dtype=float)
    np.divide(
        x_to_y_counts,
        x_core_sizes[np.newaxis, :],
        out=x_to_y,
        where=x_core_sizes[np.newaxis, :] != 0,
    )
    if not supress_progress:
        printProgressBar(
            rows // 2, rows, prefix="Progress:", suffix="Complete", length=40
        )

    # core Y against expanded X already has the public (Y, X) orientation.
    y_to_x_counts = _sketch_intersection_counts(mod_set_y, mod_set_x_neighbors)
    y_to_x = np.zeros((rows, cols), dtype=float)
    np.divide(
        y_to_x_counts,
        y_core_sizes[:, np.newaxis],
        out=y_to_x,
        where=y_core_sizes[:, np.newaxis] != 0,
    )
    symmetric_containment = np.maximum(x_to_y, y_to_x)
    containment_matrix = _identity_matrix_from_containment(
        symmetric_containment, identity, k
    )

    if not supress_progress:
        printProgressBar(
            rows, rows, prefix="Progress:", suffix="Completed", length=40
        )  # show completed progress bar
        print("\n")
    return containment_matrix


# Function used to find matching color palette to those available in const.py
def findElementsWithPrefix(lst, prefix):
    matching_elements = []
    for element in lst:
        if element.startswith(prefix):
            matching_elements.append(element)
    return matching_elements


def getInteractiveColor(palette_name, palette_orientation):
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


def getMatchingColors(color_name):
    available_colors = [
        element
        for sublist in [DIVERGING_PALETTES, QUALITATIVE_PALETTES, SEQUENTIAL_PALETTES]
        for element in sublist
    ]
    matching_elements = findElementsWithPrefix(available_colors, color_name)
    return matching_elements[-1]


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


def verifyModimizers(sparsity, l):
    # Get the next highest power of 2, if not provided
    updated_sparsity = nextPowerOfTwo(sparsity)

    sparsity_layers = [updated_sparsity]
    while l > 0:
        if sparsity_layers[-1] == 1:
            return sparsity_layers
        elif sparsity_layers[-1] % 2 == 1:
            sparsity_layers[-1] = int(sparsity_layers[-1] + 1)
        sparsity_layers.append(int(sparsity_layers[-1] / 2))
        l -= 1

    return sparsity_layers


def nextPowerOfTwo(n):
    if n <= 0:
        return 1
    n -= 1
    n |= n >> 1
    n |= n >> 2
    n |= n >> 4
    n |= n >> 8
    n |= n >> 16
    return n + 1


def generateDictionaryFromList(lst: List[int]) -> Dict[Tuple[int, int], int]:
    result = {}
    for i in range(len(lst) - 1):
        result[(lst[i], lst[i + 1])] = i
    return result


def findValueInRange(integer: int, range_dict: dict) -> int:
    if integer > max(key[0] for key in range_dict.keys()):
        return 0
    highest_value = max(range_dict.values()) + 1
    for key, value in range_dict.items():
        if key[0] >= integer >= key[1]:
            return value
    return highest_value


def setZoomLevels(axis_length, sparsity_layers):
    zoom_levels = []
    zoom_levels.append(axis_length)
    for i in range(1, len(sparsity_layers)):
        zoom_levels.append(round(axis_length / pow(2, i)))
    return zoom_levels


def makeDifferencesEqual(x, x_prime, y, y_prime):
    difference_x = abs(x_prime - x)
    difference_y = abs(y_prime - y)

    if difference_x != difference_y:
        if difference_x < difference_y:
            x_prime += difference_y - difference_x
        else:
            y_prime += difference_x - difference_y

    return x_prime, y_prime
