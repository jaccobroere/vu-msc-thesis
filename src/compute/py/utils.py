import os
import re
import uuid


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


def rename_files_double_underscore(folder_path: str):
    for filename in os.listdir(folder_path):
        print(f"Old filename: {filename}", end=" ")
        # Check if the file name contains 'estimate', 'rmsfe', or 'y_pred' preceded by an underscore
        if re.search(r"_estimate|_rmsfe|_y_pred", filename):
            # Replace the single underscore with a double underscore
            new_filename = re.sub(r"(_estimate|_rmsfe|_y_pred)", r"_\1", filename)
            # Rename the file
            os.rename(
                os.path.join(folder_path, filename),
                os.path.join(folder_path, new_filename),
            )
        else:
            print("No match")


def walk_and_replace(path: str):
    for root, dirs, files in os.walk(path):
        for direc in dirs:
            if not is_uuid(direc):
                continue
            p = os.path.join(root, direc)
            rename_files_double_underscore(p)


def is_uuid(string: str):
    try:
        uuid.UUID(string)
        return True
    except ValueError:
        return False
