from typing import Protocol, TypeVar, runtime_checkable

T = TypeVar('T', bound='Comparable')


@runtime_checkable
class Comparable(Protocol):
    def __lt__(self, other: T) -> bool: ...
    def __le__(self, other: T) -> bool: ...
    def __gt__(self, other: T) -> bool: ...
    def __ge__(self, other: T) -> bool: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...