# Chandra aimpoint monitor

This repo has content related to monitoring and correcting for the drift of
the Chandra observing aimpoint related to the temperature-driven alignment
changes of the ACA.

The key elements in this repo are:

- `aimpoint_mon` package which provides tools to monitor the aimpoint drift
  on a daily basis and update a web page which gets reviewed by the SS&AWG.
- `fit_aimpoint_drift_*.ipynb` notebooks which are used to fit discrete jumps
  in the aimpoint following a prolonged dwell at hot temperature (typically
  normal sun).
- Scripts and notebooks that were used in the analysis of aimpoint offsets
  and development of the production tools used to implement dynamical aimpoint
  offsets.
- `update_characteristics.py` script that was used for an interim period to
  generate a modified version of the OFLS characteristics file. This was used
  to remove ACA alignment drift from the target X-ray aimpoint.
