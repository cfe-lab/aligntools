
[![codecov](https://codecov.io/gh/cfe-lab/aligntools/branch/master/graph/badge.svg)](https://codecov.io/gh/cfe-lab/aligntools)
[![types - Mypy](https://img.shields.io/badge/types-Mypy-blue.svg)](https://github.com/python/mypy)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![License - AGPL3](https://img.shields.io/badge/license-AGPLv3-blue)](https://spdx.org/licenses/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/cfe-lab/aligntools/pulls)

# Alignment tools

This Python library provides robust tools for handling CIGAR strings and related alignment formats in bioinformatics. CIGAR ("Compact Idiosyncratic Gapped Alignment Report") strings succinctly represent alignment data between sequences - such as genomic sequences - highlighting matches, mismatches, insertions, deletions, and more within aligned sequences. The functionality within this library abstracts and simplifies interactions with CIGAR strings, making it easier to manipulate, interpret, and convert these alignment specifications for analytical and visualization purposes.

## Features

- **Easy CIGAR String Parsing and Construction**: Easily parse CIGAR strings from text or create them programmatically with detailed annotations for operations like matching, inserting, or deleting.
- **Alignment Manipulation**: Slice, trim, and concatenate sequence alignments for advanced genomic analyses.
- **Coordinate Mapping**: Provides bi-directional mapping between reference and query sequences, enabling you to track and modify coordinates with precision.
- **Gap Identification**: Detect gaps, insertions, and deletions within alignments to facilitate adjustments and further analysis.
- **Multiple Sequence Alignment (MSA)**: Convert CIGAR strings into multiple sequence alignment (MSA) representations for visualization or further processing.
- **CIGAR String Normalization**: Merge consecutive identical operations and support merging non-overlapping alignments to normalize complex CIGAR strings.

## Examples

### Parsing and Manipulating a CIGAR String

You can easily parse a CIGAR string and perform various operations using the `Cigar` class:

```python
from aligntools import Cigar

# Parse a CIGAR string
cigar = Cigar.parse("10M1I5D5M")

# Enumerate operations
for operation in cigar.iterate_operations():
    print(operation)

# Output:
# M M M M M M M M M M I D D D D D M M M M M

# Calculate the lengths of the reference and query sequences
print(f"Reference Length: {cigar.ref_length}")
print(f"Query Length: {cigar.query_length}")

# Output:
# Reference Length: 20
# Query Length: 16
```

### Cutting a CIGAR String by Reference Position

CIGAR strings can be split at specific reference positions, allowing for more fine-grained control over the alignment:

```python
from aligntools import Cigar, CigarHit

# Parse a CIGAR string and create a CigarHit
cigar = Cigar.coerce("10M5I10M")
hit = CigarHit(cigar, r_st=0, r_ei=19, q_st=0, q_ei=24)

# Cut the alignment at reference position 10.5 (at midpoint between positions 10 and 11).
left, right = hit.cut_reference(10.5)

print("Left slice:", left)
print("Right slice:", right)

# Output:
# Left slice: 10M5I1M@[0,15]->[0,10]
# Right slice: 9M@[16,24]->[11,19]
```

### Trimming Query and Reference Sequences

You can trim unmatched regions from either the query or reference sequence:

```python
from aligntools import Cigar

# Parse a CIGAR string with unmatched regions
cigar = Cigar.coerce("5S10M5S")

# Trim unmatched regions from the query sequence
trimmed_cigar = cigar.rstrip_query()

print(f"Trimmed CIGAR: {trimmed_cigar}")

# Output:
# Trimmed CIGAR: 5S10M
```

### Converting CIGAR Strings to Multiple Sequence Alignments (MSA)

Convert a CIGAR string to a multiple sequence alignment (MSA) representation for better visualization of how sequences align:

```python
from aligntools import Cigar

# Parse a CIGAR string
cigar = Cigar.coerce("5M2I5M")

# Define reference and query sequences
ref_seq =   "ACGTACGTAC"
query_seq = "ACGTTACGTATG"

# Convert to MSA
ref_msa, query_msa = cigar.to_msa(ref_seq, query_seq)

print(f"Reference MSA: {ref_msa}")
print(f"Query MSA:     {query_msa}")

# Output:
# Reference MSA: ACGTA--CGTAC
# Query MSA:     ACGTTACGTATG
```

### Merging Consecutive Alignments

`aligntools` can merge two consecutive CIGAR strings into a single normalized CIGAR string:

```python
from aligntools import Cigar

# Parse two CIGAR strings
cigar1 = Cigar.coerce("5M5D")
cigar2 = Cigar.coerce("10M")

# Merge the two alignments
merged_cigar = cigar1 + cigar2

print(f"Merged CIGAR: {merged_cigar}")
# Output: Merged CIGAR: 5M5D10M
```

### Advanced Usage: Coordinate Mapping Between Reference and Query Sequences

You can manage coordinate translations between the reference and query using `CoordinateMapping`:

```python
from aligntools import Cigar

# Create a CIGAR and its coordinate mapping
cigar = Cigar.coerce("5M2I3M")
mapping = cigar.coordinate_mapping

# Translate reference and query coordinates
ref_coordinate = 3
query_coordinate = mapping.ref_to_query[ref_coordinate]

print(f"Query coordinate: {query_coordinate}")
# Output:
# Query coordinate: 3

ref_coordinate = 6
query_coordinate = mapping.ref_to_query[ref_coordinate]
print(f"Query coordinate: {query_coordinate}")
# Output:
# Query coordinate: 8
```

### Using `CigarHit` for Complex Alignment Manipulations

The `CigarHit` class is a higher-level abstraction that allows you to perform complex operations on a sequence alignment:

```python
from aligntools import Cigar, CigarHit

# Define a complex CIGAR string and create a CigarHit
cigar = Cigar.coerce("5M2D5M")
hit = CigarHit(cigar, r_st=0, r_ei=11, q_st=0, q_ei=9)

# Cut and trim the alignment
left, right = hit.cut_reference(4.5)
trimmed_hit = right.rstrip_reference()

print("Left Hit:           ", left)
print("Right (trimmed) Hit:", trimmed_hit)

# Output:
# Left Hit:            5M@[0,4]->[0,4]
# Right (trimmed) Hit: 2D5M@[5,9]->[5,11]
```

## Installation

To use `aligntools` for your projects, simply run `pip install aligntools`.

## License

This project is licensed under the AGPLv3.0 License. See the COPYING file for details.
