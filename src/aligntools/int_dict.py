
from typing import Optional, AbstractSet, Mapping, \
    Dict, Iterator, MutableSet

from abc import ABC, abstractmethod, abstractproperty


class FrozenIntDict(ABC, Mapping[int, int]):
    @abstractproperty
    def domain(self) -> AbstractSet[int]:
        ...

    @abstractproperty
    def codomain(self) -> AbstractSet[int]:
        ...

    @abstractmethod
    def left_max(self, index: int) -> Optional[int]:
        ...

    @abstractmethod
    def right_min(self, index: int) -> Optional[int]:
        ...

    @abstractmethod
    def translate(self, domain_delta: int, codomain_delta: int) \
            -> 'IntDict':
        ...


class IntDict(FrozenIntDict):
    """
    An extension of the basic Python dictionary designed for
    integer-to-integer mappings.

    The IntDict maintains not just key-value pairs (as in a normal
    dictionary) but also tracks additional sets called `domain` and
    `codomain`. These sets are supersets of the keys and values
    respectively, as they include integers that might not be used
    directly in mappings but are within the range of interest for the
    domain and codomain.
    """

    def __init__(self) -> None:
        self._dict: Dict[int, int] = {}
        self._domain: MutableSet[int] = set()   # superset of self.keys()
        self._codomain: MutableSet[int] = set()  # superset of self.values()

    def extend(self, key: Optional[int], value: Optional[int]) -> None:
        if key is not None and value is not None:
            self._dict[key] = value

        if key is not None:
            self._domain.add(key)

        if value is not None:
            self._codomain.add(value)

    def left_max(self, index: int) -> Optional[int]:
        return max((v for (k, v) in self.items() if k <= index), default=None)

    def right_min(self, index: int) -> Optional[int]:
        return min((v for (k, v) in self.items() if k >= index), default=None)

    def translate(self, domain_delta: int, codomain_delta: int) -> 'IntDict':
        """
        Generates a new IntDict by shifting the entire mapping -- keys
        and values are incremented by domain_delta and codomain_delta,
        respectively.  This shift operation preserves the inherent
        ordering and relative spacing within the mapping, effectively
        repositioning the dataset within the integer space.
        """

        ret = IntDict()

        for k, v in self.items():
            ret.extend(k + domain_delta, v + codomain_delta)

        for k in self.domain:
            ret.extend(k + domain_delta, None)

        for v in self.codomain:
            ret.extend(None, v + codomain_delta)

        return ret

    @property
    def domain(self) -> AbstractSet[int]:
        return self._domain

    @property
    def codomain(self) -> AbstractSet[int]:
        return self._codomain

    def __len__(self) -> int:
        return len(self._dict)

    def __iter__(self) -> Iterator[int]:
        return iter(self._dict)

    def __getitem__(self, key: int) -> int:
        return self._dict[key]

    def __str__(self) -> str:
        pairs = []

        for k in sorted(self.domain):
            pairs.append(f"{k}: {self._dict.get(k)}")

        values_set = set(self.values())
        for v in sorted(self.codomain):
            if v not in values_set:
                pairs.append(f"None: {v}")

        return "{ " + ", ".join(pairs) + " }"
