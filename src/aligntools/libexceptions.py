

class CigarError(Exception):
    """Base class for all cigar-related custom exceptions."""
    pass


class CoersionError(TypeError, CigarError):
    """Exception raised for errors in coercing objects to CIGAR strings."""
    pass


class ParseError(ValueError, CigarError):
    """Exception raised for errors in parsing CIGAR strings."""
    pass


class InvalidOperationError(ValueError, CigarError):
    """Exception raised for invalid operations within CIGAR strings."""
    pass


class CigarHitRangeError(ValueError, CigarError):
    """Exception raised when a CIGAR hit range does not match with its coordinates."""
    pass


class CigarConnectError(ValueError, CigarError):
    """Exception raised for errors in connecting CIGAR hits."""
    pass


class CigarCutError(IndexError, CigarError):
    """Exception raised for errors in cutting CIGAR strings at reference points."""
    pass


class EmptyCigarHitListError(ValueError, CigarError):
    """Exception raised when an empty list of CIGAR hits is given for processing."""
    pass
