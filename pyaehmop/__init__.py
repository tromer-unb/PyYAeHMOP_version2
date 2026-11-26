"""Independent YAeHMOP implementation in Python (scientific engine, never invokes C)."""
__version__ = "0.4.0"
from .input import read_input
from .engine import calculate
from .friendly import read_easy_job, calculate_easy_job
