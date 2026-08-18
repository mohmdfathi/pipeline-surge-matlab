# pipeline-surge-matlab

This project studies valve-closure transients in a representative industrial
cooling-water pipeline. A steady operating point is first obtained from the
intersection of the pump and hydraulic system curves. The
resulting flow and head distribution initialize a reservoir-pipe-valve
Method-of-Characteristics model for further studies.

Several linear valve-closing times are compared using pressure histories, 
spatial pressure envelopes, the Joukowsky relation, vapour-pressure margin,
and thin-wall hoop stress. The case data are synthetic but physically 
plausible and are intended as a transparent engineering study rather than 
the design of a particular installation.

![Valve pressure history for each closing time](valve_pressure_history.jpg)

## Contents

- `analysis_valve_timing.mlx` — live script: steady operating point, MOC
  transient setup, valve-closing-time screening (pressure history, pressure
  envelopes, vapor-pressure margin, hoop stress)
- `moc_solver.m` — Method-of-Characteristics solver for the reservoir-pipe-valve
  transient
- `LICENSE` — MIT
