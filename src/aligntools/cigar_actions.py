
from enum import IntEnum

import aligntools.exceptions as ex


class CigarActions(IntEnum):
    """
    Mapping as defined on page 8 of
      <https://samtools.github.io/hts-specs/SAMv1.pdf>.
    """

    # Alignment match (can be a sequence match or mismatch)
    MATCH = 0

    # Insertion to the reference
    INSERT = 1

    # Deletion from the reference
    DELETE = 2

    # Skipped region from the reference
    SKIPPED = 3

    # Soft clip on the read (not aligned but present in the read)
    SOFT_CLIPPED = 4

    # Hard clip on the read (not present in the read)
    HARD_CLIPPED = 5

    # Padding (silent deletion from padded reference)
    PADDING = 6

    # Sequence match
    SEQ_MATCH = 7

    # Sequence mismatch
    MISMATCH = 8

    @staticmethod
    def parse(value: str) -> 'CigarActions':
        ret = OP_MAPPING.get(value)
        if ret is None:
            raise ex.ParseError(f"Invalid action: {value!r}.")
        return ret

    def __str__(self) -> str:
        return OP_MAPPING_REV[self]

    def __repr__(self) -> str:
        return OP_REPR[self]

    def relax(self) -> 'CigarActions':
        if self == CigarActions.SEQ_MATCH:
            return CigarActions.MATCH
        elif self == CigarActions.MISMATCH:
            return CigarActions.MATCH
        else:
            return self


OP_MAPPING = {
    'M': CigarActions.MATCH,
    'I': CigarActions.INSERT,
    'D': CigarActions.DELETE,
    'N': CigarActions.SKIPPED,
    'S': CigarActions.SOFT_CLIPPED,
    'H': CigarActions.HARD_CLIPPED,
    'P': CigarActions.PADDING,
    '=': CigarActions.SEQ_MATCH,
    'X': CigarActions.MISMATCH,
}


OP_REPR = {
    CigarActions.MATCH: 'CigarActions.MATCH',
    CigarActions.INSERT: 'CigarActions.INSERT',
    CigarActions.DELETE: 'CigarActions.DELETE',
    CigarActions.SKIPPED: 'CigarActions.SKIPPED',
    CigarActions.SOFT_CLIPPED: 'CigarActions.SOFT_CLIPPED',
    CigarActions.HARD_CLIPPED: 'CigarActions.HARD_CLIPPED',
    CigarActions.PADDING: 'CigarActions.PADDING',
    CigarActions.SEQ_MATCH: 'CigarActions.SEQ_MATCH',
    CigarActions.MISMATCH: 'CigarActions.MISMATCH',
}


OP_MAPPING_REV = {v: k for (k, v) in OP_MAPPING.items()}
