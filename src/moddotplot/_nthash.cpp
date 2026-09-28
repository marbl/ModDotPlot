#ifndef Py_LIMITED_API
#define Py_LIMITED_API 0x03080000
#endif
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdint>
#include <cstring>
#include <limits>

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
