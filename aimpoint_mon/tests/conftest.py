import os

# mica (imported by update_observed_aimpoints) reads $SKA at import time. The plotting
# tests never touch SKA data, so any path works when it is not already set.
os.environ.setdefault("SKA", "/proj/sot/ska")
