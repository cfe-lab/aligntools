
[![codecov](https://codecov.io/gh/cfe-lab/aligntools/branch/master/graph/badge.svg)](https://codecov.io/gh/cfe-lab/aligntools)
[![types - Mypy](https://img.shields.io/badge/types-Mypy-blue.svg)](https://github.com/python/mypy)
[![flake8 checked](https://img.shields.io/badge/flake8-checked-blueviolet.svg)](https://github.com/PyCQA/flake8)
[![License - GPL3](https://img.shields.io/badge/license-GPLv3-blue)](https://spdx.org/licenses/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/cfe-lab/aligntools/pulls)

# Alignment tools

This Python library provides robust tools for handling CIGAR strings and related alignment formats in bioinformatics. CIGAR ("Compact Idiosyncratic Gapped Alignment Report") strings succinctly represent alignment data between sequences - such as genomic sequences - highlighting matches, mismatches, insertions, deletions, and more within aligned sequences. The functionality within this library abstracts and simplifies interactions with CIGAR strings, making it easier to manipulate, interpret, and convert these alignment specifications for analytical and visualization purposes.

## Features

- **Intuitive CIGAR String Parsing and Creation**: Easily parse CIGAR strings from text formats or construct them programmatically with detailed operations annotations.
- **Alignment Manipulation**: Offers methods for slicing, trimming, and concatenating alignments, providing flexibility in processing sequence alignments.
- **Coordinate Mapping**: Leverage bidirectional mapping between reference and query sequences for detailed analysis.
- **Alignment Operations Enumeration**: Enumerate through detailed alignment operations with reference and query indexes included, facilitating deeper analysis of alignments.
- **Gap Identification**: Identify gaps, insertions, and deletions within alignments for precise manipulations or adjustments.
- **Multiple Sequence Alignment (MSA) Conversion**: Convert CIGAR string alignments into MSA representations, displaying aligned sequences with gaps.
- **CIGAR String Normalization and Merging**: Normalize CIGAR strings by merging consecutive identical operations and support merging non-overlapping but consecutive CIGAR alignments.

## Quick Start

To get started with this CIGAR string library, you can parse a CIGAR string, and execute various operations on it.

```python
from aligntools import Cigar, CigarHit, CigarActions

# Parsing a CIGAR string
cigar = Cigar.coerce("10M1I5D5M")

# Enumerate operations
for operation in cigar.iterate_operations():
    print(operation)

# Calculate reference and query lengths
print("Reference Length:", cigar.ref_length)
print("Query Length:", cigar.query_length)

# Convert to Multiple Sequence Alignment (MSA)
ref_seq = "ACGTACGTAC"
query_seq = "ACGT---CGTAC"
ref_msa, query_msa = cigar.to_msa(ref_seq, query_seq)
print("Reference MSA:", ref_msa)
print("Query MSA:", query_msa)
```

## Advanced Usage

For more sophisticated analyses, you can manipulate CIGAR strings and alignments extensively using the provided classes and methods.

### CoordinateMapping

This class allows detailed tracking and conversions between reference and query sequence coordinates, facilitating complex alignment manipulations.

```python
# Create a coordinate mapping from a CIGAR string
mapping = cigar.coordinate_mapping
```

### CigarHit

For applications that require detailed manipulation of sequence alignments—including cutting, trimming, and extending alignments—the `CigarHit` class provides numerous utilities.

```python
# Create a CigarHit and slice it at a specific reference position
hit = CigarHit(cigar, 0, 9, 0, 9)  # Define the hit range
left, right = hit.cut_reference(4.5)  # Cut at position 4.5 in reference space

print("Left Slice:", left)
print("Right Slice:", right)
```

## Contributing

Contributions to improve this library or extend its capabilities are highly welcome. Please feel free to fork the repository, make your changes, and submit a pull request.

## License

This project is licensed under the GPLv3.0 License. See the COPYING file for details.
