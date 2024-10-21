"""
Module for handling CIGAR hits.
"""

from math import ceil, floor
from typing import Tuple, Iterable, Optional, List
from dataclasses import dataclass
from functools import cached_property, reduce
from fractions import Fraction
import re

from aligntools.coordinate_mapping import CoordinateMapping
from aligntools.cigar_actions import CigarActions
from aligntools.cigar import Cigar
import aligntools.exceptions as ex


def intervals_overlap(x: Tuple[int, int], y: Tuple[int, int]) -> bool:
    """ Check if two intervals [x0, x1] and [y0, y1] overlap. """
    return x[0] <= y[1] and x[1] >= y[0]


parse_expr = re.compile(r'(?P<cigar>.+)@'
                        r'\[(?P<q_st>\d+),(?P<q_ei>\d+)\]->'
                        r'\[(?P<r_st>\d+),(?P<r_ei>\d+)\]')


@dataclass(frozen=True)
class CigarHit:
    """
    This class provides an abstraction over the complex details
    involved in working with sequence alignments expressed as CIGAR
    strings. It implements operations for alignment handling that are
    conceptually straightforward but challenging to implement ad-hoc.

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
            raise ex.CigarHitRangeError(
                f"CIGAR string maps {self.cigar.ref_length}"
                " reference positions, but CIGAR hit range"
                f" is {self.ref_length}.")

        if self.query_length != self.cigar.query_length:
            raise ex.CigarHitRangeError(
                f"CIGAR string maps {self.cigar.query_length}"
                f" query positions,"
                f" but CIGAR hit range is {self.query_length}.")

    @property
    def r_en(self):
        return self.r_ei + 1

    @property
    def q_en(self):
        return self.q_ei + 1

    @property
    def ref_length(self):
        return self.r_en - self.r_st

    @property
    def query_length(self):
        return self.q_en - self.q_st

    @staticmethod
    def from_default_alignment(r_st: int, r_ei: int, q_st: int, q_ei: int) \
            -> 'CigarHit':
        """
        A convenience method that creates a CigarHit instance
        representing a default alignment, where there are only
        deletions in the reference sequence and only insertions in the
        query.
        """

        ref_length = r_ei - r_st + 1
        query_length = q_ei - q_st + 1
        cigar = Cigar.coerce([(ref_length, CigarActions.DELETE),
                              (query_length, CigarActions.INSERT)])

        return CigarHit(cigar, r_st=r_st, r_ei=r_ei, q_st=q_st, q_ei=q_ei)

    def overlaps_in_query(self, other: 'CigarHit') -> bool:
        """
        Determines whether this CigarHit overlaps with another in
        terms of query coordinates.
        """

        return intervals_overlap((self.q_st, self.q_ei),
                                 (other.q_st, other.q_ei))

    def overlaps_in_reference(self, other: 'CigarHit') -> bool:
        """
        Determines whether this CigarHit overlaps with another in
        terms of reference coordinates.
        """

        return intervals_overlap((self.r_st, self.r_ei),
                                 (other.r_st, other.r_ei))

    def touches_in_query(self, other: 'CigarHit') -> bool:
        """
        Checks if the end of this CigarHit is immediately adjacent to
        the start of another one.
        """

        return self.q_ei + 1 == other.q_st

    def touches_in_reference(self, other: 'CigarHit') -> bool:
        """
        Checks if the end of this CigarHit is immediately adjacent to
        the start of another one.
        """

        return self.r_ei + 1 == other.r_st

    def _gaps(self, is_deletions: bool) -> Iterable['CigarHit']:
        last_query_index = self.q_st
        last_ref_index = self.r_st
        gap_start: Optional[int] = None
        op_to_ref = {v: k for k, v
                     in self.coordinate_mapping.ref_to_op.items()}
        op_to_query = {v: k for k, v
                       in self.coordinate_mapping.query_to_op.items()}
        present = op_to_ref if is_deletions else op_to_query
        missing = op_to_query if is_deletions else op_to_ref

        for op_index in sorted(self.coordinate_mapping.query_to_op.codomain) \
                + [None]:
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
                    yield CigarHit(cigar,
                                   q_st=q_st, q_ei=q_ei,
                                   r_st=r_st, r_ei=r_ei)
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
        The addition is simply a concatenation of two Cigar strings,
        and adjustment of hit coordinates.
        """

        if not (self.touches_in_query(other) and
                self.touches_in_reference(other)):
            raise ex.CigarConnectError("Cannot combine CIGAR hits that"
                                       " do not touch in both reference"
                                       " and query coordinates.")

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
            raise ex.CigarConnectError(
                "Cannot combine overlapping CIGAR hits.")

        filler = CigarHit.from_default_alignment(self.r_ei + 1, other.r_st - 1,
                                                 self.q_ei + 1, other.q_st - 1)
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
        Splits this CigarHit into two non-overlapping parts using a
        fractional cut point in the reference space.
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
        right = self._slice(left.r_ei + 1, left.q_ei + 1, ceil(op_fcut_point),
                            self.cigar.op_length)

        return left, right

    def lstrip_query(self) -> 'CigarHit':
        """
        Return a copy of the CigarHit with leading (unmatched)
        query elements removed.
        """

        cigar = self.cigar.lstrip_query()
        return CigarHit(
            cigar,
            r_st=self.r_ei - cigar.ref_length + 1, r_ei=self.r_ei,
            q_st=self.q_ei - cigar.query_length + 1, q_ei=self.q_ei)

    def rstrip_query(self) -> 'CigarHit':
        """
        Return a copy of the CigarHit with trailing (unmatched)
        query elements removed.
        """

        cigar = self.cigar.rstrip_query()
        return CigarHit(
            cigar,
            r_st=self.r_st, r_ei=self.r_st + cigar.ref_length - 1,
            q_st=self.q_st, q_ei=self.q_st + cigar.query_length - 1)

    def lstrip_reference(self) -> 'CigarHit':
        """
        Return a copy of the CigarHit with leading (unmatched)
        reference elements removed.
        """

        cigar = self.cigar.lstrip_reference()
        return CigarHit(
            cigar,
            r_st=self.r_ei - cigar.ref_length + 1, r_ei=self.r_ei,
            q_st=self.q_ei - cigar.query_length + 1, q_ei=self.q_ei)

    def rstrip_reference(self) -> 'CigarHit':
        """
        Return a copy of the CigarHit with trailing (unmatched)
        reference elements removed.
        """

        cigar = self.cigar.rstrip_reference()
        return CigarHit(
            cigar,
            r_st=self.r_st, r_ei=self.r_st + cigar.ref_length - 1,
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
        Constructs a multiple sequence alignment (MSA) representation
        for this CigarHit, using the original reference and query
        sequences. It aligns the sequences according to the CIGAR
        operations, introducing gaps ('-') as necessary to reflect
        insertions or deletions.
        """

        return self.cigar.to_msa(reference_seq[(self.r_st - 1):self.r_ei],
                                 query_seq[(self.q_st - 1):self.q_ei])

    def translate(self, reference_delta: int, query_delta: int) -> 'CigarHit':
        return CigarHit(cigar=self.cigar,
                        r_st=self.r_st + reference_delta,
                        r_ei=self.r_ei + reference_delta,
                        q_st=self.q_st + query_delta,
                        q_ei=self.q_ei + query_delta)

    @staticmethod
    def parse(string: str) -> 'CigarHit':
        """
        Parses a string representation of a CigarHit
        and returns a CigarHit object.

        This method is inverse of CigarHit.__str__.

        :param hit_str: A string representation of a CigarHit.
        :return: CigarHit object parsed from the input string.
        :raises ParseError: If the string cannot be parsed into a CigarHit.
        """

        # Regular expression to match the structure of a serialized CigarHit
        match = parse_expr.match(string)

        if not match:
            raise ex.ParseError(f"Invalid CigarHit string format: {string!r}.")

        # Extracting components from the matched regex groups
        cigar_str = match.group('cigar')
        q_st = int(match.group('q_st'))
        q_ei = int(match.group('q_ei'))
        r_st = int(match.group('r_st'))
        r_ei = int(match.group('r_ei'))

        # Validating that start indices
        # are less than or equal to end indices.
        if q_st > q_ei + 1:
            raise ex.ParseError(
                f"Query start index ({q_st}) "
                f"greater than end index ({q_ei} + 1) in: {string!r}.")

        if r_st > r_ei + 1:
            raise ex.ParseError(
                f"Reference start index ({r_st}) "
                f"greater than end index ({r_ei} + 1) in: {string!r}.")

        cigar: Cigar = Cigar.coerce(cigar_str)
        return CigarHit(cigar, r_st, r_ei, q_st, q_ei)

    def __repr__(self):
        return 'CigarHit(%r, r_st=%r, r_ei=%r, q_st=%r, q_ei=%r)' \
            % (self.cigar, self.r_st, self.r_ei, self.q_st, self.q_ei)

    def __str__(self):
        return '%s@[%d,%d]->[%d,%d]' \
            % (str(self.cigar), self.q_st, self.q_ei, self.r_st, self.r_ei)


def connect_cigar_hits(cigar_hits: List[CigarHit]) -> List[CigarHit]:
    """
    This function exists to deal with the fact that mappy does not
    always connect big gaps, and returns surrounding parts as two
    separate alignment hits.

    For those cases we simply connect all the parts that do not overlap.

    Order of cigar_hits matters because we ignore alignments
    that overlap (in query space) with previously found alignments.
    """

    if len(cigar_hits) == 0:
        raise ex.EmptyCigarHitListError(
            "Expected a non-empty list of cigar hits.")

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

    # Collect all intervals back together,
    # connecting them with CigarActions.DELETE.
    return [reduce(CigarHit.connect, group)
            for group in sorted_groups]
