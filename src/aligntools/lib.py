"""
Module for handling CIGAR strings and related alignment formats.
"""

from enum import IntEnum
from math import ceil, floor
import re
from typing import Tuple, Iterable, Optional, Set, List, Union, Dict
from dataclasses import dataclass
from functools import cached_property, reduce
from fractions import Fraction

import aligntools.libexceptions as ex


class IntDict(Dict[int, int]):
    """
    An extension of the basic Python dictionary designed for integer-to-integer mappings.

    The IntDict maintains not just key-value pairs (as in a normal dictionary) but also
    tracks additional sets called `domain` and `codomain`. These sets are supersets
    of the keys and values respectively, as they include integers that might not be used
    directly in mappings but are within the range of interest for the domain and codomain.
    """

    def __init__(self) -> None:
        super().__init__()
        self.domain: Set[int] = set()   # superset of self.keys()
        self.codomain: Set[int] = set()  # superset of self.values()

    def extend(self, key: Optional[int], value: Optional[int]) -> None:
        if key is not None and value is not None:
            self[key] = value

        if key is not None:
            self.domain.add(key)

        if value is not None:
            self.codomain.add(value)

    def left_max(self, index: int) -> Optional[int]:
        return max((v for (k, v) in self.items() if k <= index), default=None)

    def right_min(self, index: int) -> Optional[int]:
        return min((v for (k, v) in self.items() if k >= index), default=None)

    def translate(self, domain_delta: int, codomain_delta: int) -> 'IntDict':
        """
        Generates a new IntDict by shifting the entire mapping -- keys and values
        are incremented by domain_delta and codomain_delta, respectively.
        This shift operation preserves the inherent ordering and relative spacing within the mapping,
        effectively repositioning the dataset within the integer space.
        """

        ret = IntDict()

        for k, v in self.items():
            ret.extend(k + domain_delta, v + codomain_delta)

        for k in self.domain:
            ret.extend(k + domain_delta, None)

        for v in self.codomain:
            ret.extend(None, v + codomain_delta)

        return ret


class CoordinateMapping:
    """
    Manages bidirectional mappings between reference and query coordinates, as well as operation indices.

    The mapping enables conversion from reference to query coordinates and vice versa. It also manages the
    association of these coordinates with their respective operations in the alignment process.
    """

    def __init__(self) -> None:
        self.ref_to_query = IntDict()
        self.query_to_ref = IntDict()
        self.ref_to_op = IntDict()
        self.query_to_op = IntDict()

    def extend(self,
               ref_index: Optional[int],
               query_index: Optional[int],
               op_index: int) -> None:

        self.ref_to_query.extend(ref_index, query_index)
        self.query_to_ref.extend(query_index, ref_index)
        self.ref_to_op.extend(ref_index, op_index)
        self.query_to_op.extend(query_index, op_index)

    def translate(self, reference_delta: int, query_delta: int) -> 'CoordinateMapping':
        """
        Generate a new CoordinateMapping with shifted coordinate spaces.

        This method creates a new mapping where each original coordinate in
        the reference and query sequences is shifted. This allows for adapting
        the CoordinateMapping to account for changes or offsets in sequence positions,
        such as when sequences are trimmed or extended.
        """

        ret = CoordinateMapping()

        ret.ref_to_query = self.ref_to_query.translate(reference_delta, query_delta)
        ret.query_to_ref = self.query_to_ref.translate(query_delta, reference_delta)
        ret.ref_to_op = self.ref_to_op.translate(reference_delta, 0)
        ret.query_to_op = self.query_to_op.translate(query_delta, 0)

        return ret

    def __eq__(self, other):
        return (self.ref_to_op, self.query_to_op) \
            == (other.ref_to_op, other.query_to_op)

    def __repr__(self):
        return f'CoordinateMapping({self.ref_to_op},{self.query_to_op})'


# Mapping as defined in https://samtools.github.io/hts-specs/SAMv1.pdf, page 8
CigarActions = IntEnum(
    'CigarActions',
    'MATCH INSERT DELETE SKIPPED SOFT_CLIPPED HARD_CLIPPED PADDING SEQ_MATCH MISMATCH',
    start=0)


class Cigar:
    """
    Represents an alignment between a query sequence and a reference sequence using the
    Compact Idiosyncratic Gapped Alignment Report (CIGAR) string format.

    A CIGAR string is a sequence of operation codes ('M', 'I', 'D', etc.) each preceded by
    the number of bases or residues to which the operation applies.

    The class abstracts a CIGAR string as a sequence of discrete operations for convenient
    manipulation (as seen in self.iterate_operations()), while retaining the compact
    form for storage and return purposes (seen in self.__str__()).

    Instances of this class should be created by calling the `Cigar.coerce` method.
    Examples:
        Cigar.coerce("10M1I5M1D")
        Cigar.coerce([(10, CigarActions.MATCH), (1, CigarActions.INSERT), ...])
        Cigar.coerce(existing_cigar_object)

    CIGAR strings are defined in the SAM specification (https://samtools.github.io/hts-specs/SAMv1.pdf).
    """

    def __init__(self, data: Iterable[Tuple[int, CigarActions]]) -> None:
        self._data: List[Tuple[int, CigarActions]] = list(Cigar.normalize(data))

    @staticmethod
    def coerce(obj: Union['Cigar', str, Iterable[Tuple[int, CigarActions]]]) -> 'Cigar':
        if isinstance(obj, Cigar):
            return obj

        if isinstance(obj, str):
            return Cigar.parse(obj)

        if isinstance(obj, list) or isinstance(obj, tuple):
            return Cigar(obj)

        raise ex.CoersionError(f"Cannot coerce {obj!r} to CIGAR string.")

    def iterate_operations(self) -> Iterable[CigarActions]:
        """
        Yields each operation in the CIGAR sequence as a `CigarActions` enum.
        The resulting sequence is a decoded version of the initial run-length encoded sequence.
        """

        for num, operation in self._data:
            for _ in range(num):
                yield operation

    def iterate_operations_with_pointers(self) -> Iterable[Tuple[CigarActions, Optional[int], Optional[int]]]:
        """
        Iterates over the operations while tracking the reference and
        query sequence positions affected by each operation.

        Example:
            For a Cigar instance representing "1M1I1M", this method would yield:
            (CigarActions.MATCH, 0, 0), (CigarActions.INSERT, None, 1), (CigarActions.MATCH, 1, 2)

        :return: Tuple of type (CigarActions, reference_pointer, query_pointer) for each operation in the
            CIGAR sequence. Pointers can be None if the operation does not map to a sequence
            position (e.g., insertions, deletions).
        """

        ref_pointer = 0
        query_pointer = 0

        for operation in self.iterate_operations():
            if operation in (CigarActions.MATCH, CigarActions.SEQ_MATCH, CigarActions.MISMATCH):
                yield operation, ref_pointer, query_pointer
                query_pointer += 1
                ref_pointer += 1

            elif operation in (CigarActions.INSERT, CigarActions.SOFT_CLIPPED):
                yield operation, None, query_pointer
                query_pointer += 1

            elif operation in (CigarActions.DELETE, CigarActions.SKIPPED):
                yield operation, ref_pointer, None
                ref_pointer += 1

            else:
                yield operation, None, None

    def slice_operations(self, start_inclusive: int, end_noninclusive: int) -> 'Cigar':
        """
        Creates a new Cigar object by slicing the current one from start_inclusive to
        end_noninclusive. Note that slicing is done at the level of individual operations,
        not at the level of counts within operations.

        Example:
            Given a Cigar instance representing "10M5D5M", slicing from 2 to 11 would result in a new
            Cigar object representing "8M1D".
        """

        return Cigar([(1, op) for op in self.iterate_operations()]
                     [start_inclusive:end_noninclusive])

    def lstrip_query(self) -> 'Cigar':
        """ Return a copy of the Cigar with leading (unmatched) query elements removed. """

        min_r = min(self.coordinate_mapping.ref_to_query.keys(), default=None)
        if min_r is None:
            min_op = float("inf")
        else:
            min_op = self.coordinate_mapping.ref_to_op.get(min_r, float("inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if query_pointer is None or i >= min_op]
        return Cigar.coerce(ops)

    def rstrip_query(self) -> 'Cigar':
        """ Return a copy of the Cigar with trailing (unmatched) query elements removed. """

        max_r = max(self.coordinate_mapping.ref_to_query.keys(), default=None)
        if max_r is None:
            max_op = float("-inf")
        else:
            max_op = self.coordinate_mapping.ref_to_op.get(max_r, float("-inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if query_pointer is None or i <= max_op]
        return Cigar.coerce(ops)

    def lstrip_reference(self) -> 'Cigar':
        """ Return a copy of the Cigar with leading (unmatched) reference elements removed. """

        min_q = min(self.coordinate_mapping.query_to_ref.keys(), default=None)
        if min_q is None:
            min_op = float("inf")
        else:
            min_op = self.coordinate_mapping.query_to_op.get(min_q, float("inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if ref_pointer is None or i >= min_op]
        return Cigar.coerce(ops)

    def rstrip_reference(self) -> 'Cigar':
        """ Return a copy of the Cigar with trailing (unmatched) reference elements removed. """

        max_q = max(self.coordinate_mapping.query_to_ref.keys(), default=None)
        if max_q is None:
            max_op = float("-inf")
        else:
            max_op = self.coordinate_mapping.query_to_op.get(max_q, float("-inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if ref_pointer is None or i <= max_op]
        return Cigar.coerce(ops)

    @cached_property
    def coordinate_mapping(self) -> CoordinateMapping:
        """
        Convert this CIGAR string to coordinate mapping representing
        a reference-to-query and query-to-reference coordinate mappings.

        :return: Lists of integers representing the mappings of coordinates from the reference
                sequence to the query sequence, and back.
        """

        mapping = CoordinateMapping()

        for op_pointer, (operation, ref_pointer, query_pointer) in enumerate(self.iterate_operations_with_pointers()):
            mapping.extend(ref_pointer,
                           query_pointer,
                           op_pointer)

        return mapping

    def to_msa(self, reference_seq: str, query_seq: str) -> Tuple[str, str]:
        """
        Constructs a multiple sequence alignment (MSA) representation for this Cigar, using the original reference
        and query sequences. It aligns the sequences according to the CIGAR operations, introducing gaps ('-')
        as necessary to reflect insertions or deletions.
        """

        reference_msa = ''
        query_msa = ''

        for operation, ref_pointer, query_pointer in self.iterate_operations_with_pointers():
            if ref_pointer is None and query_pointer is None:
                continue

            try:
                if ref_pointer is not None:
                    reference_msa += reference_seq[ref_pointer]
                else:
                    reference_msa += '-'

                if query_pointer is not None:
                    query_msa += query_seq[query_pointer]
                else:
                    query_msa += '-'

            except IndexError:
                raise ex.ParseError("CIGAR string corresponds to a larger match than either reference or query.")

        return reference_msa, query_msa

    @staticmethod
    def from_msa(reference: str, query: str) -> 'Cigar':
        """
        Converts a Multiple Sequence Alignment (MSA) of a reference and a query sequence
        into a CIGAR object. Alignments are expected to be strings of equal length where
        gaps are represented by '-' characters.

        :param reference: The reference sequence in MSA format.
        :param query: The query sequence in MSA format.
        :return: A Cigar object representing the alignment.
        """

        if len(reference) != len(query):
            raise ex.ParseError("Reference and query sequences must be of the same length.")

        operations = []
        curr_op = ''
        count = 0

        for ref_base, query_base in zip(reference, query):
            if ref_base == '-' and query_base == '-':
                # This scenario should not happen in a correct MSA
                continue
            elif ref_base == '-':
                op = 'I'  # Insertion in reference
            elif query_base == '-':
                op = 'D'  # Deletion from reference
            elif ref_base == query_base:
                op = 'M'  # Match/Mismatch (Exact match)
            else:
                op = 'M'  # Match/Mismatch (Mismatch)

            if op == curr_op:
                count += 1
            else:
                if curr_op:  # Not the first operation
                    operations.append((count, Cigar.parse_operation(curr_op)))
                curr_op = op
                count = 1

        # Don't forget to add the last operation
        if curr_op:
            operations.append((count, Cigar.parse_operation(curr_op)))

        return Cigar(operations)

    @cached_property
    def op_length(self):
        return sum(1 for x in self.iterate_operations())

    @cached_property
    def query_length(self):
        return max((query_pointer + 1 if query_pointer is not None else 0 for (_, _, query_pointer)
                    in self.iterate_operations_with_pointers()),
                   default=0)

    @cached_property
    def ref_length(self):
        return max((ref_pointer + 1 if ref_pointer is not None else 0 for (_, ref_pointer, _)
                    in self.iterate_operations_with_pointers()),
                   default=0)

        #                                 #
        #  Boring boilerplate code below  #
        #                                 #

    OP_MAPPING = {
        'M': CigarActions.MATCH,         # Alignment match (can be a sequence match or mismatch)
        'I': CigarActions.INSERT,        # Insertion to the reference
        'D': CigarActions.DELETE,        # Deletion from the reference
        'N': CigarActions.SKIPPED,       # Skipped region from the reference
        'S': CigarActions.SOFT_CLIPPED,  # Soft clip on the read (ignored region, not aligned but present in the read)
        'H': CigarActions.HARD_CLIPPED,  # Hard clip on the read (ignored region, not present in the read)
        'P': CigarActions.PADDING,       # Padding (silent deletion from padded reference, not applicable for our case)
        '=': CigarActions.SEQ_MATCH,     # Sequence match
        'X': CigarActions.MISMATCH,      # Sequence mismatch
    }

    @staticmethod
    def parse_operation(operation: str) -> CigarActions:
        if operation in Cigar.OP_MAPPING:
            return Cigar.OP_MAPPING[operation]
        else:
            raise ex.InvalidOperationError(f"Unexpected CIGAR action: {operation!r}.")

    @staticmethod
    def operation_to_str(op: CigarActions) -> str:
        return [k for (k, v) in Cigar.OP_MAPPING.items() if v == op][0]

    @staticmethod
    def parse(string: str) -> 'Cigar':
        """
        Parses a CIGAR string into a Cigar object.

        :param string: A CIGAR string with the format '(\\d+[MIDNSHPX=])+', where each operation code
                       is preceded by a number indicating how many times the operation should be applied.
        """

        data = []
        while string:
            match = re.match(r'([0-9]+)([^0-9])', string)
            if match:
                num, operation = match.groups()
                data.append((int(num), Cigar.parse_operation(operation)))
                string = string[match.end():]
            else:
                raise ex.ParseError(f"Invalid CIGAR string. Invalid part: {string[:20]!r}.")

        return Cigar(data)

    @staticmethod
    def normalize(cigar_lst: Iterable[Tuple[int, CigarActions]]) -> Iterable[Tuple[int, CigarActions]]:
        """
        Goes through the list appending operations to the CIGAR sequence,
        checking for type correctness and performing normalization
        by merging consecutive identical operations.
        """

        last_item = None

        for item in cigar_lst:
            # Type checking
            if not isinstance(item, list) and not isinstance(item, tuple):
                raise ex.InvalidOperationError(f"Invalid CIGAR list: {item!r} is not a tuple.")
            if len(item) != 2:
                raise ex.InvalidOperationError(f"Invalid CIGAR list: {item!r} has a bad length.")

            num, operation = item
            if isinstance(operation, int):
                operation = CigarActions(operation)
            if not isinstance(num, int) or not isinstance(operation, CigarActions):
                raise ex.InvalidOperationError(f"Invalid CIGAR list: {item!r} is not a number/operation tuple.")
            if num < 0:
                raise ex.InvalidOperationError("Invalid CIGAR list: number of operations is negative.")

            # Normalization
            if num == 0:
                continue

            if last_item:
                last_num, last_operation = last_item
                if operation == last_operation:
                    last_item = (last_num + num, operation)
                    continue

            if last_item:
                yield last_item[0], last_item[1]
            last_item = item

        if last_item:
            yield last_item[0], last_item[1]

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Cigar) and self._data == other._data

    def __add__(self, other: 'Cigar') -> 'Cigar':
        return Cigar(self._data + other._data)

    def __repr__(self):
        return f'Cigar({str(self)!r})'

    def __str__(self):
        """ Inverse of `Cigar.parse` """
        return ''.join('{}{}'.format(num, Cigar.operation_to_str(op)) for num, op in self._data)


def intervals_overlap(x: Tuple[int, int], y: Tuple[int, int]) -> bool:
    """ Check if two intervals [x0, x1] and [y0, y1] overlap. """
    return x[0] <= y[1] and x[1] >= y[0]


@dataclass(frozen=True)
class CigarHit:
    """
    This class provides an abstraction over the complex details involved in working with sequence alignments
    expressed as CIGAR strings. It implements operations for alignment handling that are conceptually
    straightforward but challenging to implement ad-hoc.

    The main tasks handled by this class are:
    - Precisely dividing an alignment into two contiguous segments
      at any given reference position (`cut_reference()`),
    - Removing portions of the query sequence that do not align with
      the reference sequence from either end
      while preserving the alignment context (`lstrip*()` and `rstrip*()`),
    - Enumerating gaps in the alignment (`gaps()`).
    """

    cigar: Cigar
    r_st: int
    r_ei: int  # inclusive
    q_st: int
    q_ei: int  # inclusive

    def __post_init__(self):
        if self.ref_length != self.cigar.ref_length:
            raise ex.CigarHitRangeError(f"CIGAR string maps {self.cigar.ref_length}"
                                        f" reference positions, but CIGAR hit range is {self.ref_length}.")

        if self.query_length != self.cigar.query_length:
            raise ex.CigarHitRangeError(f"CIGAR string maps {self.cigar.query_length}"
                                        f" query positions, but CIGAR hit range is {self.query_length}.")

    @property
    def ref_length(self):
        return self.r_ei + 1 - self.r_st

    @property
    def query_length(self):
        return self.q_ei + 1 - self.q_st

    @staticmethod
    def from_default_alignment(r_st: int, r_ei: int, q_st: int, q_ei: int) -> 'CigarHit':
        """
        A convenience method that creates a CigarHit instance representing a default alignment,
        where there are only deletions in the reference sequence and only insertions in the query.
        """

        ref_length = r_ei - r_st + 1
        query_length = q_ei - q_st + 1
        cigar = Cigar.coerce([(ref_length, CigarActions.DELETE),
                              (query_length, CigarActions.INSERT)])

        return CigarHit(cigar, r_st=r_st, r_ei=r_ei, q_st=q_st, q_ei=q_ei)

    def overlaps_in_query(self, other: 'CigarHit') -> bool:
        """
        Determines whether this CigarHit overlaps with another in terms of query coordinates.
        """

        return intervals_overlap((self.q_st, self.q_ei), (other.q_st, other.q_ei))

    def overlaps_in_reference(self, other: 'CigarHit') -> bool:
        """
        Determines whether this CigarHit overlaps with another in terms of reference coordinates.
        """

        return intervals_overlap((self.r_st, self.r_ei), (other.r_st, other.r_ei))

    def touches_in_query(self, other: 'CigarHit') -> bool:
        """
        Checks if the end of this CigarHit is immediately adjacent to the start of another one.
        """

        return self.q_ei + 1 == other.q_st

    def touches_in_reference(self, other: 'CigarHit') -> bool:
        """
        Checks if the end of this CigarHit is immediately adjacent to the start of another one.
        """

        return self.r_ei + 1 == other.r_st

    def _gaps(self, is_deletions: bool) -> Iterable['CigarHit']:
        last_query_index = self.q_st
        last_ref_index = self.r_st
        gap_start: Optional[int] = None
        op_to_ref = {v: k for k, v in self.coordinate_mapping.ref_to_op.items()}
        op_to_query = {v: k for k, v in self.coordinate_mapping.query_to_op.items()}
        present = op_to_ref if is_deletions else op_to_query
        missing = op_to_query if is_deletions else op_to_ref

        for op_index in sorted(self.coordinate_mapping.query_to_op.codomain) + [None]:
            if op_index in present and \
               op_index not in missing:
                if gap_start is None:
                    gap_start = op_index
            else:
                if gap_start is not None:
                    if op_index is None:
                        op_index = -1

                    cigar = self.cigar.slice_operations(gap_start, op_index)
                    if is_deletions:
                        q_st = last_query_index
                        r_st = present[gap_start]
                    else:
                        q_st = present[gap_start]
                        r_st = last_ref_index
                    q_ei = q_st + cigar.query_length - 1
                    r_ei = r_st + cigar.ref_length - 1
                    yield CigarHit(cigar, q_st=q_st, q_ei=q_ei, r_st=r_st, r_ei=r_ei)
                    gap_start = None

            if op_index in op_to_query:
                last_query_index = op_to_query[op_index]
            if op_index in op_to_ref:
                last_ref_index = op_to_ref[op_index]

    def deletions(self) -> Iterable['CigarHit']:
        return self._gaps(is_deletions=True)

    def insertions(self) -> Iterable['CigarHit']:
        return self._gaps(is_deletions=False)

    def __add__(self, other: 'CigarHit') -> 'CigarHit':
        """
        Only adds CigarHits that are touching.
        The addition is simply a concatenation of two Cigar strings, and adjustment of hit coordinates.
        """

        if not (self.touches_in_query(other) and self.touches_in_reference(other)):
            raise ex.CigarConnectError("Cannot combine CIGAR hits that do not touch "
                                       "in both reference and query coordinates.")

        return CigarHit(cigar=self.cigar + other.cigar,
                        r_st=self.r_st,
                        r_ei=other.r_ei,
                        q_st=self.q_st,
                        q_ei=other.q_ei)

    def connect(self, other: 'CigarHit') -> 'CigarHit':
        """
        Inserts deletions/insertions between self and other,
        then adjusts boundaries appropriately.
        """

        if self.overlaps_in_query(other) or self.overlaps_in_reference(other):
            raise ex.CigarConnectError("Cannot combine overlapping CIGAR hits.")

        filler = CigarHit.from_default_alignment(self.r_ei + 1, other.r_st - 1, self.q_ei + 1, other.q_st - 1)
        return self + filler + other

    @property
    def epsilon(self):
        return Fraction(1, self.cigar.op_length * 3 + 1)

    def _ref_cut_to_op_cut(self, cut_point: Fraction) -> Fraction:
        mapping = self.coordinate_mapping

        left_op_cut_point = mapping.ref_to_op.left_max(floor(cut_point))
        right_op_cut_point = mapping.ref_to_op.right_min(ceil(cut_point))

        if left_op_cut_point is None:
            left_op_cut_point = -1
        if right_op_cut_point is None:
            right_op_cut_point = self.cigar.op_length

        def lerp(start: int, end: int, t: Fraction) -> Fraction:
            return (1 - t) * start + t * end

        op_cut_point = lerp(left_op_cut_point, right_op_cut_point,
                            cut_point - floor(cut_point))

        if float(op_cut_point).is_integer():
            # Disambiguate to the right.
            op_cut_point += self.epsilon

        return op_cut_point

    def _slice(self, r_st: int, q_st: int, o_st: int, o_ei: int) -> 'CigarHit':
        cigar = self.cigar.slice_operations(o_st, o_ei + 1)
        r_ei = r_st + cigar.ref_length - 1
        q_ei = q_st + cigar.query_length - 1

        return CigarHit(cigar=cigar,
                        r_st=r_st,
                        r_ei=r_ei,
                        q_st=q_st,
                        q_ei=q_ei,
                        )

    def cut_reference(self, cut_point: float) -> Tuple['CigarHit', 'CigarHit']:
        """
        Splits this CigarHit into two non-overlapping parts using a fractional cut point in the reference space.
        Resulting parts of CigarHits are touching at cut point.
        The two parts do not share any elements, and no element is "lost".
        """

        fcut_point: Fraction = Fraction(cut_point)
        if fcut_point.denominator == 1:
            raise ex.CigarCutError("Cut accepts fractions, not integers.")

        if self.ref_length == 0 or \
           not (self.r_st - 1 < fcut_point < self.r_ei + 1):
            raise ex.CigarCutError("Cut point out of reference bounds.")

        op_fcut_point = self._ref_cut_to_op_cut(fcut_point)
        left = self._slice(self.r_st, self.q_st, 0, floor(op_fcut_point))
        right = self._slice(left.r_ei + 1, left.q_ei + 1, ceil(op_fcut_point), self.cigar.op_length)

        return left, right

    def lstrip_query(self) -> 'CigarHit':
        """ Return a copy of the CigarHit with leading (unmatched) query elements removed. """

        cigar = self.cigar.lstrip_query()
        return CigarHit(cigar, r_st=self.r_ei - cigar.ref_length + 1, r_ei=self.r_ei,
                        q_st=self.q_ei - cigar.query_length + 1, q_ei=self.q_ei)

    def rstrip_query(self) -> 'CigarHit':
        """ Return a copy of the CigarHit with trailing (unmatched) query elements removed. """

        cigar = self.cigar.rstrip_query()
        return CigarHit(cigar, r_st=self.r_st, r_ei=self.r_st + cigar.ref_length - 1,
                        q_st=self.q_st, q_ei=self.q_st + cigar.query_length - 1)

    def lstrip_reference(self) -> 'CigarHit':
        """ Return a copy of the CigarHit with leading (unmatched) reference elements removed. """

        cigar = self.cigar.lstrip_reference()
        return CigarHit(cigar, r_st=self.r_ei - cigar.ref_length + 1, r_ei=self.r_ei,
                        q_st=self.q_ei - cigar.query_length + 1, q_ei=self.q_ei)

    def rstrip_reference(self) -> 'CigarHit':
        """ Return a copy of the CigarHit with trailing (unmatched) reference elements removed. """

        cigar = self.cigar.rstrip_reference()
        return CigarHit(cigar, r_st=self.r_st, r_ei=self.r_st + cigar.ref_length - 1,
                        q_st=self.q_st, q_ei=self.q_st + cigar.query_length - 1)

    @cached_property
    def coordinate_mapping(self) -> CoordinateMapping:
        """
        Convert this alignment to coordinate mapping representing
        a reference-to-query and query-to-reference coordinate mappings.
        """

        return self.cigar.coordinate_mapping.translate(self.r_st, self.q_st)

    def to_msa(self, reference_seq: str, query_seq: str) -> Tuple[str, str]:
        """
        Constructs a multiple sequence alignment (MSA) representation for this CigarHit, using the original reference
        and query sequences. It aligns the sequences according to the CIGAR operations, introducing gaps ('-')
        as necessary to reflect insertions or deletions.
        """

        return self.cigar.to_msa(reference_seq[self.r_st:], query_seq[self.q_st:])

    def translate(self, reference_delta: int, query_delta: int) -> 'CigarHit':
        return CigarHit(cigar=self.cigar,
                        r_st=self.r_st + reference_delta,
                        r_ei=self.r_ei + reference_delta,
                        q_st=self.q_st + query_delta,
                        q_ei=self.q_ei + query_delta)

    def __repr__(self):
        return 'CigarHit(%r, r_st=%r, r_ei=%r, q_st=%r, q_ei=%r)' \
            % (self.cigar, self.r_st, self.r_ei, self.q_st, self.q_ei)

    def __str__(self):
        return '%s@[%d,%d]->[%d,%d]' \
            % (str(self.cigar), self.q_st, self.q_ei, self.r_st, self.r_ei)


def connect_cigar_hits(cigar_hits: List[CigarHit]) -> List[CigarHit]:
    """
    This function exists to deal with the fact that mappy does not always
    connect big gaps, and returns surrounding parts as two separate alignment hits.

    For those cases we simply connect all the parts that do not overlap.

    Order of cigar_hits matters because we ignore alignments
    that overlap (in query space) with previously found alignments.
    """

    if len(cigar_hits) == 0:
        raise ex.EmptyCigarHitListError("Expected a non-empty list of cigar hits.")

    accumulator: List[CigarHit] = []

    # Collect non-overlapping parts.
    # Earlier matches have priority over ones that come after.
    for hit in cigar_hits:
        if any(earlier.overlaps_in_query(hit) for earlier in accumulator):
            continue

        accumulator.append(hit)

    # Sort by interval start positions.
    sorted_parts = sorted(accumulator, key=lambda p: p.r_st)

    # Segregate independent matches.
    sorted_groups: List[List[CigarHit]] = []

    def find_group(phit: CigarHit) -> None:
        for group in sorted_groups:
            if phit.q_st > group[-1].q_ei and \
               all(not phit.overlaps_in_reference(other) for other in group):
                group.append(phit)
                return

        sorted_groups.append([phit])

    for hit in sorted_parts:
        find_group(hit)

    # Collect all intervals back together, connecting them with CigarActions.DELETE.
    return [reduce(lambda x, y: x.connect(y), group) for group in sorted_groups]
