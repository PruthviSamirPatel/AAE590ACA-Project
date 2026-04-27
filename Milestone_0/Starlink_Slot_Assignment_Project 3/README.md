# Starlink Slot Assignment + Reduced PMP + Lyapunov Control

This is a self-contained MATLAB project built around the functions already used in the original project.  
It creates a Starlink-like shell, defines 5 launched satellites and 5 open slots in a RAAN/mean-anomaly map, assigns satellites to slots, and then runs both:

1. a reduced-order path-constrained PMP planner, and  
2. a Lyapunov SMA/RAAN-tracking controller.

## Run

```matlab
RUN_STARLINK_SLOT_PROJECT
```

## Modeling assumptions

The reduced model keeps only the variables requested for this project:

\[
x = [a,\Omega,M]^T
\]

with fixed circularity and fixed inclination. The control is low-thrust tangential acceleration, so

\[
\dot a = \frac{2}{n(a)}u_T,\qquad |u_T|\le u_{\max},
\]

while the RAAN evolves naturally from the secular \(J_2\) drift:

\[
\dot\Omega =
-\frac{3}{2}J_2 n(a)\left(\frac{R_E}{a}\right)^2\cos i .
\]

The PMP solution uses the minimum-time structure for this reduced problem:

- stay at the low parking/drift SMA to accumulate RAAN precession,
- then raise SMA at maximum tangential thrust.

That is the path-constrained PMP structure for the case where lower RAAN slots are reached by exploiting faster \(J_2\) precession at lower altitude.

The Lyapunov controller tracks the same reference structure using a saturated feedback law in semi-major axis.

## Main output

The script prints:

- Starlink-like scenario parameters,
- initial satellite RAAN/mean-anomaly map,
- open slot RAAN/mean-anomaly map,
- cost matrix,
- selected satellite-to-slot assignment,
- PMP timing, coast duration, raise duration, delta-V, RAAN error, and mean-anomaly residual,
- Lyapunov timing, delta-V, RAAN error, and mean-anomaly residual.

It also saves:

```text
results/starlink_assignment_results.mat
results/assignment_summary.csv
```

## Edit scenario

Edit:

```matlab
Starlink_Functions/create_starlink_scenario.m
```

to change:

- launch epoch,
- parking altitude,
- final shell altitude,
- inclination,
- thrust acceleration,
- slot RAAN offsets,
- slot mean anomalies.
