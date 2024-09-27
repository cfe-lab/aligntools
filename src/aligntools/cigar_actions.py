
from enum import IntEnum


class CigarActions(IntEnum):
    """
    Mapping as defined on page 8 of
      <https://samtools.github.io/hts-specs/SAMv1.pdf>.
    """

    MATCH = 0
    INSERT = 1
    DELETE = 2
    SKIPPED = 3
    SOFT_CLIPPED = 4
    HARD_CLIPPED = 5
    PADDING = 6
    SEQ_MATCH = 7
    MISMATCH = 8
