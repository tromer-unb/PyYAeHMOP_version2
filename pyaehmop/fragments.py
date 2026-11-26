"""Atom range syntax in FMO fragments: '1-5,8,10-12'."""
import re


def parse_integer_ranges(raw):
    result=[]
    for token in raw.replace(" ", "").split(","):
        if not token:
            continue
        match=re.fullmatch(r"(\d+)(?:-(\d+))?",token)
        if match is None:
            raise ValueError(f"Invalid FMO fragment: {token!r}")
        low=int(match.group(1))
        high=int(match.group(2)) if match.group(2) else low
        if low<1 or high<low:
            raise ValueError(f"Invalid atom range: {token!r}")
        result.extend(range(low,high+1))
    if not result or len(result)!=len(set(result)):
        raise ValueError("Empty fragment or duplicated atoms")
    return result
