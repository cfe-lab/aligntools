"""
Module for handling CIGAR strings.
"""

import re
from typing import Tuple, Iterable, Optional, Union
from functools import cached_property

from aligntools.coordinate_mapping import CoordinateMapping
from aligntools.cigar_actions import CigarActions
import aligntools.exceptions as ex


class Cigar:
    """
    Represents an alignment between a query sequence and a reference
    sequence using the Compact Idiosyncratic Gapped Alignment Report
    (CIGAR) string format.

    A CIGAR string is a sequence of operation codes ('M', 'I', 'D', etc.)
    each preceded by the number of bases or residues to which
    the operation applies.

    The class abstracts a CIGAR string as a sequence of discrete
    operations for convenient manipulation (as seen in
    self.iterate_operations()), while retaining the compact form for
    storage and return purposes (seen in self.__str__()).

    Instances of this class should be created
    by calling the `Cigar.coerce` method.
    Examples:
        Cigar.coerce("10M1I5M1D")
        Cigar.coerce([(10, CigarActions.MATCH), (1, CigarActions.INSERT), ...])
        Cigar.coerce(existing_cigar_object)

    CIGAR strings are defined in the SAM specification
       (https://samtools.github.io/hts-specs/SAMv1.pdf).
    """

    def __init__(self, data: Iterable[Tuple[int, CigarActions]]) -> None:
        self._data: Tuple[Tuple[int, CigarActions], ...] \
            = tuple(Cigar.normalize(data))

    @staticmethod
    def coerce(obj: Union['Cigar', str, Iterable[Tuple[int, CigarActions]]]) \
            -> 'Cigar':
        if isinstance(obj, Cigar):
            return obj

        if isinstance(obj, str):
            return Cigar.parse(obj)

        if isinstance(obj, list) or isinstance(obj, tuple):
            return Cigar(obj)

        raise ex.CoersionError(f"Cannot coerce {obj!r} to CIGAR string.")

    def iterate_operations(self) -> Iterable[CigarActions]:
        """
        Yields each operation in the CIGAR sequence as a
        `CigarActions` enum.  The resulting sequence is a decoded
        version of the initial run-length encoded sequence.
        """

        for num, operation in self._data:
            for _ in range(num):
                yield operation

    def iterate_operations_with_pointers(self) \
            -> Iterable[Tuple[CigarActions, Optional[int], Optional[int]]]:
        """
        Iterates over the operations while tracking the reference and
        query sequence positions affected by each operation.

        Example:
            For a Cigar instance representing "1M1I1M",
            this method would yield:
            (CigarActions.MATCH, 0, 0),
            (CigarActions.INSERT, None, 1),
            (CigarActions.MATCH, 1, 2)

        :return: Tuple of type (CigarActions, reference_pointer, query_pointer)
            for each operation in the CIGAR sequence. Pointers can be None
            if the operation does not map to a sequence
            position (e.g., insertions, deletions).
        """

        ref_pointer = 0
        query_pointer = 0

        for operation in self.iterate_operations():
            if operation in (CigarActions.MATCH,
                             CigarActions.SEQ_MATCH,
                             CigarActions.MISMATCH):
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

    def slice_operations(self, start_inclusive: int, end_noninclusive: int) \
            -> 'Cigar':
        """
        Creates a new Cigar object by slicing the current one from
        start_inclusive to end_noninclusive. Note that slicing is done
        at the level of individual operations, not at the level of
        counts within operations.

        Example:
            Given a Cigar instance representing "10M5D5M",
            slicing from 2 to 11 would result in a new
            Cigar object representing "8M1D".
        """

        return Cigar([(1, op) for op in self.iterate_operations()]
                     [start_inclusive:end_noninclusive])

    def lstrip_query(self) -> 'Cigar':
        """
        Return a copy of the Cigar with leading (unmatched)
        query elements removed.
        """

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
        """
        Return a copy of the Cigar with trailing (unmatched)
        query elements removed.
        """

        max_r = max(self.coordinate_mapping.ref_to_query.keys(), default=None)
        if max_r is None:
            max_op = float("-inf")
        else:
            max_op = self.coordinate_mapping.ref_to_op.get(
                max_r, float("-inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if query_pointer is None or i <= max_op]
        return Cigar.coerce(ops)

    def lstrip_reference(self) -> 'Cigar':
        """
        Return a copy of the Cigar with leading (unmatched)
        reference elements removed.
        """

        min_q = min(self.coordinate_mapping.query_to_ref.keys(), default=None)
        if min_q is None:
            min_op = float("inf")
        else:
            min_op = self.coordinate_mapping.query_to_op.get(
                min_q, float("inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if ref_pointer is None or i >= min_op]
        return Cigar.coerce(ops)

    def rstrip_reference(self) -> 'Cigar':
        """
        Return a copy of the Cigar with trailing (unmatched)
        reference elements removed.
        """

        max_q = max(self.coordinate_mapping.query_to_ref.keys(), default=None)
        if max_q is None:
            max_op = float("-inf")
        else:
            max_op = self.coordinate_mapping.query_to_op.get(
                max_q, float("-inf"))

        ops = [(1, op) for i, (op, ref_pointer, query_pointer)
               in enumerate(self.iterate_operations_with_pointers())
               if ref_pointer is None or i <= max_op]
        return Cigar.coerce(ops)

    @cached_property
    def coordinate_mapping(self) -> CoordinateMapping:
        """
        Convert this CIGAR string to coordinate mapping representing
        a reference-to-query and query-to-reference coordinate mappings.

        :return: Lists of integers representing the mappings of coordinates
                 from the reference sequence to the query sequence, and back.
        """

        mapping = CoordinateMapping()

        for op_pointer, (operation, ref_pointer, query_pointer) \
                in enumerate(self.iterate_operations_with_pointers()):
            mapping.extend(ref_pointer,
                           query_pointer,
                           op_pointer)

        return mapping

    def to_msa(self, reference_seq: str, query_seq: str) -> Tuple[str, str]:
        """
        Constructs a multiple sequence alignment (MSA) representation
        for this Cigar, using the original reference and query
        sequences. It aligns the sequences according to the CIGAR
        operations, introducing gaps ('-') as necessary to reflect
        insertions or deletions.
        """

        reference_msa = ''
        query_msa = ''

        for operation, ref_pointer, query_pointer \
                in self.iterate_operations_with_pointers():
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
                raise \
                    ex.MSALengthError("CIGAR string corresponds to a larger"
                                      " match than either reference or query.")

        return reference_msa, query_msa

    @staticmethod
    def from_msa(reference: str, query: str) -> 'Cigar':
        """
        Converts a Multiple Sequence Alignment (MSA) of a reference
        and a query sequence into a CIGAR object. Alignments are
        expected to be strings of equal length where gaps are
        represented by '-' characters.

        :param reference: The reference sequence in MSA format.
        :param query: The query sequence in MSA format.
        :return: A Cigar object representing the alignment.
        """

        if len(reference) != len(query):
            raise ex.ParseError("Reference and query sequences"
                                " must be of the same length.")

        operations = []

        for ref_base, query_base in zip(reference, query):
            if ref_base == '-' and query_base == '-':
                # This scenario should not happen in a correct MSA
                continue
            elif ref_base == '-':
                op = CigarActions.INSERT
            elif query_base == '-':
                op = CigarActions.DELETE
            elif ref_base == query_base:
                op = CigarActions.MATCH
            else:
                op = CigarActions.MATCH

            operations.append((1, op))

        return Cigar(operations)

    @cached_property
    def op_length(self):
        return sum(1 for x in self.iterate_operations())

    @cached_property
    def query_length(self):
        return max((query_pointer + 1
                    if query_pointer is not None else 0
                    for (_, _, query_pointer)
                    in self.iterate_operations_with_pointers()),
                   default=0)

    @cached_property
    def ref_length(self):
        return max((ref_pointer + 1
                    if ref_pointer is not None else 0
                    for (_, ref_pointer, _)
                    in self.iterate_operations_with_pointers()),
                   default=0)

        #                                 #
        #  Boring boilerplate code below  #
        #                                 #

    @staticmethod
    def parse_operation(operation: str) -> CigarActions:
        try:
            return CigarActions.parse(operation)
        except ex.ParseError:
            raise ex.InvalidOperationError(
                f"Unexpected CIGAR action: {operation!r}.")

    @staticmethod
    def parse(string: str) -> 'Cigar':
        """
        Parses a CIGAR string into a Cigar object.

        :param string: A CIGAR string with the format '(\\d+[MIDNSHPX=])+',
                       where each operation code is preceded by a number
                       indicating how many times the operation
                       should be applied.
        """

        data = []
        while string:
            match = re.match(r'([0-9]+)([^0-9])', string)
            if match:
                num, operation = match.groups()
                data.append((int(num), Cigar.parse_operation(operation)))
                string = string[match.end():]
            else:
                raise ex.ParseError("Invalid CIGAR string."
                                    f" Invalid part: {string[:20]!r}.")

        return Cigar(data)

    @staticmethod
    def normalize(cigar_lst: Iterable[Tuple[int, CigarActions]]) \
            -> Iterable[Tuple[int, CigarActions]]:
        """
        Goes through the list appending operations to the CIGAR sequence,
        checking for type correctness and performing normalization
        by merging consecutive identical operations.
        """

        last_item = None

        for item in cigar_lst:
            # Type checking
            if not isinstance(item, list) and not isinstance(item, tuple):
                raise ex.InvalidOperationError(f"Invalid CIGAR list: {item!r}"
                                               " is not a tuple.")
            if len(item) != 2:
                raise ex.InvalidOperationError(f"Invalid CIGAR list: {item!r}"
                                               " has a bad length.")

            num, operation = item
            if isinstance(operation, int):
                operation = CigarActions(operation)
            if not isinstance(num, int) or \
               not isinstance(operation, CigarActions):
                raise ex.InvalidOperationError(
                    f"Invalid CIGAR list: {item!r}"
                    " is not a number/operation tuple.")
            if num < 0:
                raise ex.InvalidOperationError(
                    "Invalid CIGAR list: number of operations is negative.")

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
        return ''.join('{}{}'.format(num, op) for num, op in self._data)
