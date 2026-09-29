#ifndef Py_LIMITED_API
#define Py_LIMITED_API 0x03080000
#endif
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>

// The rolling hash kernel below is derived from ntHash v2.4.0, commit
// c26bd4572a19de81e30d55042dbd33c1fd21d4b6. ntHash is distributed under
// the MIT license; see LICENSE.ntHash in the source distribution.

namespace {

constexpr uint64_t MASK_33 = (uint64_t{1} << 33) - 1;
constexpr uint64_t MASK_31 = (uint64_t{1} << 31) - 1;
constexpr uint64_t SEED_A = 0x3c8bfbb395c60474ULL;
constexpr uint64_t SEED_C = 0x3193c18562a02b4cULL;
constexpr uint64_t SEED_G = 0x20323ed082572324ULL;
constexpr uint64_t SEED_T = 0x295549f54be24456ULL;
constexpr uint64_t FNV_OFFSET_BASIS = 14695981039346656037ULL;
constexpr uint64_t FNV_PRIME = 1099511628211ULL;

inline char
upper_ascii(char value)
{
  return value >= 'a' && value <= 'z'
           ? static_cast<char>(value - ('a' - 'A'))
           : value;
}

inline uint64_t
seed(char value)
{
  switch (upper_ascii(value)) {
    case 'A': return SEED_A;
    case 'C': return SEED_C;
    case 'G': return SEED_G;
    case 'T':
    case 'U': return SEED_T;
    default: return 0;
  }
}

inline uint64_t
complement_seed(char value)
{
  switch (upper_ascii(value)) {
    case 'A': return SEED_T;
    case 'C': return SEED_G;
    case 'G': return SEED_C;
    case 'T':
    case 'U': return SEED_A;
    default: return 0;
  }
}

inline bool
valid(char value)
{
  return seed(value) != 0;
}

inline uint64_t
rotate_left_width(uint64_t value, unsigned amount, unsigned width)
{
  amount %= width;
  if (amount == 0) {
    return value;
  }
  const uint64_t mask = width == 33 ? MASK_33 : MASK_31;
  return ((value << amount) | (value >> (width - amount))) & mask;
}

inline uint64_t
rotate_right_width(uint64_t value, unsigned amount, unsigned width)
{
  amount %= width;
  if (amount == 0) {
    return value;
  }
  const uint64_t mask = width == 33 ? MASK_33 : MASK_31;
  return ((value >> amount) | (value << (width - amount))) & mask;
}

inline uint64_t
srol(uint64_t value, unsigned amount = 1)
{
  const uint64_t low = rotate_left_width(value & MASK_33, amount, 33);
  const uint64_t high = rotate_left_width(value >> 33, amount, 31);
  return low | (high << 33);
}

inline uint64_t
sror(uint64_t value)
{
  const uint64_t low = rotate_right_width(value & MASK_33, 1, 33);
  const uint64_t high = rotate_right_width(value >> 33, 1, 31);
  return low | (high << 33);
}

inline void
base_hashes(const char* sequence,
            unsigned k,
            uint64_t& forward,
            uint64_t& reverse)
{
  forward = 0;
  reverse = 0;
  for (unsigned i = 0; i < k; ++i) {
    forward = srol(forward) ^ seed(sequence[i]);
    reverse = srol(reverse) ^ complement_seed(sequence[k - i - 1]);
  }
}

inline uint64_t
canonical(uint64_t forward, uint64_t reverse)
{
  // Unsigned overflow supplies the modulo-2^64 operation specified by
  // ntHash2's canonical hash.
  return forward + reverse;
}

inline unsigned char
normalized_forward_char(char value)
{
  const char normalized = upper_ascii(value);
  return static_cast<unsigned char>(normalized == 'U' ? 'T' : normalized);
}

inline unsigned char
normalized_reverse_complement_char(char value)
{
  switch (upper_ascii(value)) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T':
    case 'U': return 'A';
    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'M': return 'K';
    case 'K': return 'M';
    case 'S': return 'S';
    case 'W': return 'W';
    case 'H': return 'D';
    case 'B': return 'V';
    case 'V': return 'B';
    case 'D': return 'H';
    case 'N': return 'N';
    default: return static_cast<unsigned char>(upper_ascii(value));
  }
}

inline uint64_t
fnv1a(const char* sequence, unsigned k, bool reverse_complement)
{
  uint64_t hash = FNV_OFFSET_BASIS;
  for (unsigned i = 0; i < k; ++i) {
    const unsigned char value = reverse_complement
                                  ? normalized_reverse_complement_char(
                                      sequence[k - i - 1])
                                  : normalized_forward_char(sequence[i]);
    hash ^= value;
    hash *= FNV_PRIME;
  }
  return hash;
}

void
add_interval_layer(const char* sequence,
                   Py_ssize_t start,
                   Py_ssize_t end,
                   unsigned k,
                   bool use_canonical,
                   bool include_ambiguous,
                   unsigned sparsity_power,
                   std::vector<uint64_t>& selected)
{
  if (start >= end) {
    return;
  }
  unsigned invalid_count = 0;
  for (unsigned i = 0; i < k; ++i) {
    invalid_count += !valid(sequence[start + i]);
  }
  uint64_t forward = 0;
  uint64_t reverse = 0;
  bool rolling = false;
  const uint64_t mask = sparsity_power == 0
                          ? 0
                          : (uint64_t{ 1 } << sparsity_power) - 1;

  for (Py_ssize_t position = start; position < end; ++position) {
    if (position > start) {
      const bool previous_window_was_valid = invalid_count == 0;
      invalid_count -= !valid(sequence[position - 1]);
      invalid_count += !valid(sequence[position + k - 1]);
      if (invalid_count == 0 && previous_window_was_valid) {
        const char outgoing = sequence[position - 1];
        const char incoming = sequence[position + k - 1];
        forward =
          srol(forward) ^ seed(incoming) ^ srol(seed(outgoing), k);
        reverse ^= srol(complement_seed(incoming), k);
        reverse ^= complement_seed(outgoing);
        reverse = sror(reverse);
        rolling = true;
      } else if (invalid_count != 0) {
        rolling = false;
      }
    }

    const bool ambiguous = invalid_count != 0;
    if (ambiguous && !include_ambiguous) {
      continue;
    }
    uint64_t value = 0;
    if (!ambiguous) {
      if (!rolling) {
        base_hashes(sequence + position, k, forward, reverse);
        rolling = true;
      }
      value = use_canonical ? canonical(forward, reverse) : forward;
    } else {
      const uint64_t fallback_forward = fnv1a(sequence + position, k, false);
      value = use_canonical
                ? canonical(fallback_forward,
                            fnv1a(sequence + position, k, true))
                : fallback_forward;
    }
    if ((value & mask) == 0) {
      selected.push_back(value);
    }
  }
}

// The ordinary sparsity layer is populated during the single chromosome pass.
// Only a genuinely deficient interval is re-hashed at successively denser
// layers. Human windows overwhelmingly take the one-pass path, while the rare
// fallback remains bit-for-bit equivalent and bounded to that interval.
struct WindowSketch
{
  Py_ssize_t start;
  Py_ssize_t end;
  std::size_t output_index;
  unsigned sparsity_power;
  std::size_t minimum_size;
  std::vector<uint64_t> top_layer;

  WindowSketch(Py_ssize_t interval_start,
               Py_ssize_t interval_end,
               std::size_t index,
               unsigned power,
               std::size_t minimum)
    : start(interval_start)
    , end(interval_end)
    , output_index(index)
    , sparsity_power(power)
    , minimum_size(minimum)
  {
  }

  void add(uint64_t value)
  {
    const uint64_t top_mask = sparsity_power == 0
                                ? 0
                                : (uint64_t{ 1 } << sparsity_power) - 1;
    if ((value & top_mask) == 0) {
      top_layer.push_back(value);
    }
  }

  std::vector<uint64_t> finish(const char* sequence,
                               unsigned k,
                               bool use_canonical,
                               bool include_ambiguous)
  {
    std::sort(top_layer.begin(), top_layer.end());
    top_layer.erase(
      std::unique(top_layer.begin(), top_layer.end()), top_layer.end());
    unsigned selected_power = sparsity_power;
    while (top_layer.size() < minimum_size && selected_power > 0) {
      --selected_power;
      add_interval_layer(sequence,
                         start,
                         end,
                         k,
                         use_canonical,
                         include_ambiguous,
                         selected_power,
                         top_layer);
      std::sort(top_layer.begin(), top_layer.end());
      top_layer.erase(
        std::unique(top_layer.begin(), top_layer.end()), top_layer.end());
    }
    return std::move(top_layer);
  }
};

inline bool
is_power_of_two(Py_ssize_t value)
{
  return value > 0 &&
         (static_cast<uint64_t>(value) &
          (static_cast<uint64_t>(value) - 1)) == 0;
}

inline unsigned
power_of_two_exponent(Py_ssize_t value)
{
  unsigned power = 0;
  while (value > 1) {
    value >>= 1;
    ++power;
  }
  return power;
}

PyObject*
sketch_kmers(PyObject*, PyObject* args)
{
  const char* sequence = nullptr;
  Py_ssize_t sequence_length = 0;
  Py_ssize_t requested_k = 0;
  int use_canonical = 0;
  PyObject* intervals_object = nullptr;
  Py_ssize_t sparsity = 0;
  Py_ssize_t requested_minimum_size = 0;
  int include_ambiguous = 0;
  if (!PyArg_ParseTuple(args,
                        "s#npOnnp:sketch_kmers",
                        &sequence,
                        &sequence_length,
                        &requested_k,
                        &use_canonical,
                        &intervals_object,
                        &sparsity,
                        &requested_minimum_size,
                        &include_ambiguous)) {
    return nullptr;
  }
  if (requested_k < 1 ||
      requested_k > std::numeric_limits<uint16_t>::max()) {
    PyErr_SetString(PyExc_ValueError, "k must be between 1 and 65535");
    return nullptr;
  }
  if (!is_power_of_two(sparsity)) {
    PyErr_SetString(PyExc_ValueError, "sparsity must be a positive power of two");
    return nullptr;
  }

  const Py_ssize_t count = sequence_length >= requested_k
                             ? sequence_length - requested_k + 1
                             : 0;
  const Py_ssize_t interval_count = PySequence_Size(intervals_object);
  if (interval_count < 0) {
    return nullptr;
  }
  const unsigned sparsity_power = power_of_two_exponent(sparsity);
  const std::size_t minimum_size = requested_minimum_size > 0
                                     ? static_cast<std::size_t>(
                                         requested_minimum_size)
                                     : 0;

  std::vector<WindowSketch> windows;
  std::vector<std::vector<uint64_t>> output;
  std::vector<std::size_t> schedule;
  try {
    windows.reserve(static_cast<std::size_t>(interval_count));
    for (Py_ssize_t index = 0; index < interval_count; ++index) {
      PyObject* interval = PySequence_GetItem(intervals_object, index);
      if (interval == nullptr) {
        return nullptr;
      }
      Py_ssize_t start = 0;
      Py_ssize_t end = 0;
      const int parsed = PyArg_ParseTuple(interval, "nn", &start, &end);
      Py_DECREF(interval);
      if (!parsed) {
        return nullptr;
      }
      if (start < 0 || end < start || end > count) {
        PyErr_SetString(
          PyExc_ValueError,
          "interval bounds must satisfy 0 <= start <= end <= k-mer count");
        return nullptr;
      }
      windows.emplace_back(start,
                           end,
                           static_cast<std::size_t>(index),
                           sparsity_power,
                           minimum_size);
    }

    output.resize(static_cast<std::size_t>(interval_count));
    schedule.resize(static_cast<std::size_t>(interval_count));
    std::iota(schedule.begin(), schedule.end(), std::size_t{ 0 });
    std::sort(schedule.begin(), schedule.end(), [&](std::size_t left,
                                                     std::size_t right) {
      if (windows[left].start != windows[right].start) {
        return windows[left].start < windows[right].start;
      }
      return windows[left].end < windows[right].end;
    });
  } catch (const std::bad_alloc&) {
    PyErr_NoMemory();
    return nullptr;
  } catch (const std::exception& error) {
    PyErr_SetString(PyExc_RuntimeError, error.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "failed to initialize native sketches");
    return nullptr;
  }

  std::exception_ptr computation_error;
  Py_BEGIN_ALLOW_THREADS
  try {
    std::vector<std::size_t> active;
    active.reserve(8);
    std::size_t next_window = 0;
    const unsigned k = static_cast<unsigned>(requested_k);
    unsigned invalid_count = 0;
    for (unsigned i = 0; i < k && count > 0; ++i) {
      invalid_count += !valid(sequence[i]);
    }
    uint64_t forward = 0;
    uint64_t reverse = 0;
    bool rolling = false;

    for (Py_ssize_t position = 0; position < count; ++position) {
      while (next_window < schedule.size() &&
             windows[schedule[next_window]].start <= position) {
        const std::size_t index = schedule[next_window++];
        auto& window = windows[index];
        if (window.start == window.end) {
          output[window.output_index] = window.finish(
            sequence, k, use_canonical != 0, include_ambiguous != 0);
        } else {
          active.push_back(index);
        }
      }
      auto kept_end = std::remove_if(
        active.begin(), active.end(), [&](std::size_t index) {
          auto& window = windows[index];
          if (window.end > position) {
            return false;
          }
          output[window.output_index] = window.finish(
            sequence, k, use_canonical != 0, include_ambiguous != 0);
          return true;
        });
      active.erase(kept_end, active.end());

      if (position > 0) {
        const bool previous_window_was_valid = invalid_count == 0;
        invalid_count -= !valid(sequence[position - 1]);
        invalid_count += !valid(sequence[position + requested_k - 1]);
        if (invalid_count == 0 && previous_window_was_valid) {
          const char outgoing = sequence[position - 1];
          const char incoming = sequence[position + requested_k - 1];
          forward =
            srol(forward) ^ seed(incoming) ^ srol(seed(outgoing), k);
          reverse ^= srol(complement_seed(incoming), k);
          reverse ^= complement_seed(outgoing);
          reverse = sror(reverse);
          rolling = true;
        } else if (invalid_count != 0) {
          rolling = false;
        }
      }

      const bool ambiguous = invalid_count != 0;
      if (ambiguous && !include_ambiguous) {
        continue;
      }
      uint64_t value = 0;
      if (!ambiguous) {
        if (!rolling) {
          base_hashes(sequence + position, k, forward, reverse);
          rolling = true;
        }
        value = use_canonical ? canonical(forward, reverse) : forward;
      } else {
        const uint64_t fallback_forward =
          fnv1a(sequence + position, k, false);
        value = use_canonical
                  ? canonical(fallback_forward,
                              fnv1a(sequence + position, k, true))
                  : fallback_forward;
      }
      for (const std::size_t index : active) {
        windows[index].add(value);
      }
    }

    for (const std::size_t index : active) {
      auto& window = windows[index];
      output[window.output_index] = window.finish(
        sequence, k, use_canonical != 0, include_ambiguous != 0);
    }
    while (next_window < schedule.size()) {
      auto& window = windows[schedule[next_window++]];
      output[window.output_index] = window.finish(
        sequence, k, use_canonical != 0, include_ambiguous != 0);
    }
  } catch (...) {
    computation_error = std::current_exception();
  }
  Py_END_ALLOW_THREADS

  if (computation_error != nullptr) {
    try {
      std::rethrow_exception(computation_error);
    } catch (const std::bad_alloc&) {
      PyErr_NoMemory();
    } catch (const std::exception& error) {
      PyErr_SetString(PyExc_RuntimeError, error.what());
    } catch (...) {
      PyErr_SetString(PyExc_RuntimeError, "native sketch construction failed");
    }
    return nullptr;
  }

  PyObject* result = PyList_New(interval_count);
  if (result == nullptr) {
    return nullptr;
  }
  for (Py_ssize_t index = 0; index < interval_count; ++index) {
    const auto& sketch = output[static_cast<std::size_t>(index)];
    if (sketch.size() > static_cast<std::size_t>(
                          std::numeric_limits<Py_ssize_t>::max()) /
                          sizeof(uint64_t)) {
      Py_DECREF(result);
      PyErr_SetString(PyExc_OverflowError, "sketch output is too large");
      return nullptr;
    }
    PyObject* packed = PyBytes_FromStringAndSize(
      reinterpret_cast<const char*>(sketch.data()),
      static_cast<Py_ssize_t>(sketch.size() * sizeof(uint64_t)));
    if (packed == nullptr) {
      Py_DECREF(result);
      return nullptr;
    }
    PyList_SetItem(result, index, packed);
  }
  return result;
}

struct SketchInput
{
  const char* values;
  std::size_t length;
  std::size_t row;
  bool right;
  PyObject* owner;
};

struct OwnedSketchInputs
{
  std::vector<SketchInput> values;

  OwnedSketchInputs() = default;
  OwnedSketchInputs(const OwnedSketchInputs&) = delete;
  OwnedSketchInputs& operator=(const OwnedSketchInputs&) = delete;

  ~OwnedSketchInputs()
  {
    // Destruction happens only while the calling thread owns the GIL: native
    // computation catches exceptions before Py_END_ALLOW_THREADS. Keeping one
    // reference per buffer makes the raw pointers safe even when callers pass
    // a custom sequence that creates bytes objects on demand, or mutate an
    // input list from another Python thread while the merge is running.
    for (const auto& input : values) {
      Py_DECREF(input.owner);
    }
  }
};

inline uint64_t
sketch_value_at(const SketchInput& input, std::size_t offset)
{
  uint64_t value = 0;
  std::memcpy(&value, input.values + offset * sizeof(uint64_t), sizeof(value));
  return value;
}

struct SketchCursor
{
  uint64_t value;
  std::size_t input;
  std::size_t offset;
};

struct CursorGreater
{
  bool operator()(const SketchCursor& left, const SketchCursor& right) const
  {
    return left.value > right.value;
  }
};

bool
append_sketch_inputs(PyObject* sketches,
                     Py_ssize_t count,
                     bool right,
                     OwnedSketchInputs& inputs)
{
  for (Py_ssize_t row = 0; row < count; ++row) {
    PyObject* sketch = PySequence_GetItem(sketches, row);
    if (sketch == nullptr) {
      return false;
    }
    char* buffer = nullptr;
    Py_ssize_t buffer_length = 0;
    const int status =
      PyBytes_AsStringAndSize(sketch, &buffer, &buffer_length);
    if (status < 0) {
      Py_DECREF(sketch);
      return false;
    }
    if (buffer_length < 0 ||
        buffer_length % static_cast<Py_ssize_t>(sizeof(uint64_t)) != 0) {
      Py_DECREF(sketch);
      PyErr_SetString(
        PyExc_ValueError, "sketch buffers must contain packed uint64 values");
      return false;
    }
    try {
      inputs.values.push_back(
        { buffer,
          static_cast<std::size_t>(buffer_length) / sizeof(uint64_t),
          static_cast<std::size_t>(row),
          right,
          sketch });
    } catch (const std::bad_alloc&) {
      Py_DECREF(sketch);
      PyErr_NoMemory();
      return false;
    } catch (const std::exception& error) {
      Py_DECREF(sketch);
      PyErr_SetString(PyExc_RuntimeError, error.what());
      return false;
    } catch (...) {
      Py_DECREF(sketch);
      PyErr_SetString(PyExc_RuntimeError, "failed to retain sketch input");
      return false;
    }
  }
  return true;
}

PyObject*
intersection_counts(PyObject*, PyObject* args)
{
  PyObject* left_object = nullptr;
  PyObject* right_object = nullptr;
  if (!PyArg_ParseTuple(
        args, "OO:intersection_counts", &left_object, &right_object)) {
    return nullptr;
  }
  const Py_ssize_t rows = PySequence_Size(left_object);
  if (rows < 0) {
    return nullptr;
  }
  const Py_ssize_t columns = PySequence_Size(right_object);
  if (columns < 0) {
    return nullptr;
  }
  if (columns != 0 &&
      rows > std::numeric_limits<Py_ssize_t>::max() / columns) {
    PyErr_SetString(PyExc_OverflowError, "intersection matrix is too large");
    return nullptr;
  }
  const Py_ssize_t cell_count = rows * columns;
  if (cell_count > std::numeric_limits<Py_ssize_t>::max() /
                     static_cast<Py_ssize_t>(sizeof(int32_t))) {
    PyErr_SetString(PyExc_OverflowError, "intersection matrix is too large");
    return nullptr;
  }
  if (rows > std::numeric_limits<Py_ssize_t>::max() - columns) {
    PyErr_SetString(PyExc_OverflowError, "too many sketch inputs");
    return nullptr;
  }

  OwnedSketchInputs inputs;
  try {
    inputs.values.reserve(static_cast<std::size_t>(rows + columns));
  } catch (const std::bad_alloc&) {
    PyErr_NoMemory();
    return nullptr;
  } catch (const std::exception& error) {
    PyErr_SetString(PyExc_RuntimeError, error.what());
    return nullptr;
  } catch (...) {
    PyErr_SetString(PyExc_RuntimeError, "failed to allocate sketch inputs");
    return nullptr;
  }
  // Use the dimensions captured above rather than asking arbitrary Python
  // sequences for their lengths again. A stateful __len__ must not be able to
  // produce row indices outside the matrix allocated from the first result.
  if (!append_sketch_inputs(left_object, rows, false, inputs) ||
      !append_sketch_inputs(right_object, columns, true, inputs)) {
    return nullptr;
  }

  std::vector<int32_t> counts;
  std::exception_ptr computation_error;
  Py_BEGIN_ALLOW_THREADS
  try {
    counts.assign(static_cast<std::size_t>(cell_count), int32_t{ 0 });
    std::priority_queue<SketchCursor,
                        std::vector<SketchCursor>,
                        CursorGreater>
      queue;
    for (std::size_t index = 0; index < inputs.values.size(); ++index) {
      if (inputs.values[index].length != 0) {
        queue.push({ sketch_value_at(inputs.values[index], 0), index, 0 });
      }
    }

    std::vector<std::size_t> left_rows;
    std::vector<std::size_t> right_rows;
    left_rows.reserve(static_cast<std::size_t>(rows));
    right_rows.reserve(static_cast<std::size_t>(columns));
    while (!queue.empty()) {
      const uint64_t current_hash = queue.top().value;
      left_rows.clear();
      right_rows.clear();
      do {
        const SketchCursor cursor = queue.top();
        queue.pop();
        const SketchInput& input = inputs.values[cursor.input];
        (input.right ? right_rows : left_rows).push_back(input.row);
        const std::size_t next_offset = cursor.offset + 1;
        if (next_offset < input.length) {
          queue.push(
            { sketch_value_at(input, next_offset), cursor.input, next_offset });
        }
      } while (!queue.empty() && queue.top().value == current_hash);

      for (const std::size_t left_row : left_rows) {
        const std::size_t row_offset =
          left_row * static_cast<std::size_t>(columns);
        for (const std::size_t right_row : right_rows) {
          ++counts[row_offset + right_row];
        }
      }
    }
  } catch (...) {
    computation_error = std::current_exception();
  }
  Py_END_ALLOW_THREADS

  if (computation_error != nullptr) {
    try {
      std::rethrow_exception(computation_error);
    } catch (const std::bad_alloc&) {
      PyErr_NoMemory();
    } catch (const std::exception& error) {
      PyErr_SetString(PyExc_RuntimeError, error.what());
    } catch (...) {
      PyErr_SetString(PyExc_RuntimeError,
                      "native sketch intersection failed");
    }
    return nullptr;
  }
  return PyBytes_FromStringAndSize(
    reinterpret_cast<const char*>(counts.data()),
    cell_count * static_cast<Py_ssize_t>(sizeof(int32_t)));
}

PyObject*
hash_kmers(PyObject*, PyObject* args)
{
  const char* sequence = nullptr;
  Py_ssize_t sequence_length = 0;
  Py_ssize_t requested_k = 0;
  int use_canonical = 0;
  if (!PyArg_ParseTuple(args,
                        "s#np:hash_kmers",
                        &sequence,
                        &sequence_length,
                        &requested_k,
                        &use_canonical)) {
    return nullptr;
  }

  if (requested_k < 1 ||
      requested_k > std::numeric_limits<uint16_t>::max()) {
    PyErr_SetString(PyExc_ValueError, "k must be between 1 and 65535");
    return nullptr;
  }

  const unsigned k = static_cast<unsigned>(requested_k);
  const Py_ssize_t count =
    sequence_length >= requested_k ? sequence_length - requested_k + 1 : 0;
  if (count > std::numeric_limits<Py_ssize_t>::max() /
                static_cast<Py_ssize_t>(sizeof(uint64_t))) {
    PyErr_SetString(PyExc_OverflowError, "hash output is too large");
    return nullptr;
  }

  bool has_ambiguous_base = false;
  if (count > 0) {
    Py_BEGIN_ALLOW_THREADS
    for (Py_ssize_t i = 0; i < sequence_length; ++i) {
      if (!valid(sequence[i])) {
        has_ambiguous_base = true;
        break;
      }
    }
    Py_END_ALLOW_THREADS
  }

  PyObject* hash_bytes = PyBytes_FromStringAndSize(
    nullptr, count * static_cast<Py_ssize_t>(sizeof(uint64_t)));
  if (hash_bytes == nullptr) {
    return nullptr;
  }
  PyObject* mask_bytes = PyBytes_FromStringAndSize(
    nullptr, has_ambiguous_base ? count : 0);
  if (mask_bytes == nullptr) {
    Py_DECREF(hash_bytes);
    return nullptr;
  }

  char* hash_output = PyBytes_AsString(hash_bytes);
  char* mask_output =
    has_ambiguous_base ? PyBytes_AsString(mask_bytes) : nullptr;
  if (hash_output == nullptr || (has_ambiguous_base && mask_output == nullptr)) {
    Py_DECREF(hash_bytes);
    Py_DECREF(mask_bytes);
    return nullptr;
  }

  if (count > 0) {
    Py_BEGIN_ALLOW_THREADS
    unsigned invalid_count = 0;
    for (unsigned i = 0; i < k; ++i) {
      invalid_count += !valid(sequence[i]);
    }

    uint64_t forward = 0;
    uint64_t reverse = 0;
    bool rolling = false;

    for (Py_ssize_t position = 0; position < count; ++position) {
      if (position > 0) {
        const bool previous_window_was_valid = invalid_count == 0;
        invalid_count -= !valid(sequence[position - 1]);
        invalid_count += !valid(sequence[position + requested_k - 1]);

        if (invalid_count == 0 && previous_window_was_valid) {
          const char outgoing = sequence[position - 1];
          const char incoming = sequence[position + requested_k - 1];
          forward =
            srol(forward) ^ seed(incoming) ^ srol(seed(outgoing), k);
          reverse ^= srol(complement_seed(incoming), k);
          reverse ^= complement_seed(outgoing);
          reverse = sror(reverse);
          rolling = true;
        } else if (invalid_count != 0) {
          rolling = false;
        }
      }

      const bool ambiguous = invalid_count != 0;
      uint64_t value = 0;
      if (!ambiguous) {
        if (!rolling) {
          base_hashes(sequence + position, k, forward, reverse);
          rolling = true;
        }
        value = use_canonical ? canonical(forward, reverse) : forward;
      } else {
        const uint64_t fallback_forward =
          fnv1a(sequence + position, k, false);
        value = use_canonical
                  ? canonical(fallback_forward,
                              fnv1a(sequence + position, k, true))
                  : fallback_forward;
      }

      std::memcpy(hash_output + position * sizeof(uint64_t),
                  &value,
                  sizeof(value));
      if (mask_output != nullptr) {
        mask_output[position] = static_cast<char>(ambiguous);
      }
    }
    Py_END_ALLOW_THREADS
  }

  PyObject* result = PyTuple_Pack(2, hash_bytes, mask_bytes);
  Py_DECREF(hash_bytes);
  Py_DECREF(mask_bytes);
  return result;
}

PyMethodDef methods[] = {
  { "intersection_counts",
    intersection_counts,
    METH_VARARGS,
    "intersection_counts(left, right, /)\n--\n\n"
    "Count intersections between sorted unique uint64 sketch buffers using a "
    "bounded-memory k-way merge. Return a packed row-major int32 matrix." },
  { "sketch_kmers",
    sketch_kmers,
    METH_VARARGS,
    "sketch_kmers(sequence, k, canonical, intervals, sparsity, "
    "minimum_size, include_ambiguous, /)\n--\n\n"
    "Hash a sequence and construct sorted, unique window sketches without "
    "materializing its positional hash array. intervals contains zero-based, "
    "half-open k-mer-index bounds. Return one packed native-endian uint64 "
    "buffer per interval." },
  { "hash_kmers",
    hash_kmers,
    METH_VARARGS,
    "hash_kmers(sequence, k, canonical, /)\n--\n\n"
    "Hash every positional k-mer in one native call.\n\n"
    "Return (hashes, ambiguity_mask). hashes contains packed native-endian "
    "uint64 values. ambiguity_mask is empty when every window is valid; "
    "otherwise it contains one byte per window, with 1 marking a window "
    "that contains a non-ACGTU character." },
  { nullptr, nullptr, 0, nullptr },
};

PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_nthash",
  "Stable-ABI batch bindings for the vendored ntHash2 rolling hash.",
  -1,
  methods,
};

} // namespace

PyMODINIT_FUNC
PyInit__nthash()
{
  PyObject* created_module = PyModule_Create(&module);
  if (created_module == nullptr) {
    return nullptr;
  }
  if (PyModule_AddStringConstant(
        created_module, "ALGORITHM", "ntHash_v2") < 0 ||
      PyModule_AddStringConstant(
        created_module, "UPSTREAM_VERSION", "2.4.0") < 0 ||
      PyModule_AddStringConstant(
        created_module,
        "UPSTREAM_COMMIT",
        "c26bd4572a19de81e30d55042dbd33c1fd21d4b6") < 0) {
    Py_DECREF(created_module);
    return nullptr;
  }
  return created_module;
}
