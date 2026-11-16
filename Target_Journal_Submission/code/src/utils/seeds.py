"""Deterministic seeding and hardware detection for reproducibility."""
from __future__ import annotations
import os, random
import numpy as np

def set_global_seed(seed: int) -> None:
    """Seed Python, NumPy, and PyTorch (if present) and force deterministic/thread-stable behaviour."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        torch.use_deterministic_algorithms(True, warn_only=True)
    except Exception:
        pass

def detect_device() -> str:
    """Return the best available torch device string ('cuda', 'mps', or 'cpu'); 'cpu' if torch absent."""
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda"
        if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass
    return "cpu"
