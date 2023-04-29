import re


def parse_float(s):
    """
    Parses a string and extracts the first floating-point number or scientific notation in it.

    Parameters
    ----------
    s : str
        The string to parse.

    Returns
    -------
    float or None
        The first floating-point number or scientific notation found in the string, or None if no match is found.

    Examples
    --------
    >>> parse_float('The answer is 42.')
    42.0

    >>> parse_float('1.23e-4 is a very small number.')
    0.000123

    >>> parse_float('This string contains no numbers.')
    None
    """
    pattern = r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
    match = re.search(pattern, s)
    if match:
        return float(match.group())
    else:
        return None
