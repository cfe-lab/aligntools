
from typing import Optional

from aligntools.int_dict import IntDict


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
