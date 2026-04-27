# Starlink Slot Assignment — Reduced Lyapunov Only

This project is the **reduced Lyapunov-only** version of the Starlink-like
slot assignment workflow. It removes all PMP code and uses only the reduced
state

\[
x = [a,\ \Delta\Omega,\ \Delta M]^T
\]

with circular low-thrust + \(J_2\) dynamics:

\[
\dot a = \frac{2u_T}{n(a)}, \qquad
\dot{\Delta\Omega}=\dot\Omega_{J2}(a)-\dot\Omega_{J2}(a_f),
\qquad
\dot{\Delta M}=n(a)-n(a_f).
\]

The control is scalar tangential Lyapunov feedback:

\[
V=\frac{1}{2}(a-a_f)^2,\qquad
u_T = \mathrm{sat}\{-k_a(a-a_f), \pm u_{\max}\}.
\]

RAAN targeting is handled by selecting a coast/drift time at the lower
parking shell before the Lyapunov orbit raise.

## Run

```matlab
RUN_STARLINK_REDUCED_LYAPUNOV_PROJECT
```

## Important unit note

In this reduced dynamics model, `uMax_km_s2` is an **acceleration**. The default
is restored to the Project 3 low-thrust value:

```matlab
cfg.uMax_km_s2 = 2.0e-7;   % 0.0002 m/s^2 = 0.2 mm/s^2
```

If this is changed to \(1~\mathrm{m/s^2}\), the 350 km to 550 km raise requires
only about \(112\) seconds of thrust, so the altitude history will look almost
vertical on a multi-day plot. That is a units/scale effect, not a Lyapunov
trajectory feature.


## V3 plot corrections

This version fixes the reduced-Lyapunov plots:

1. The control plot is shown in **mm/s²**, not m/s², so the low-thrust value
   \(u_{\max}=0.2\ \mathrm{mm/s^2}\) is visible for SAT01 and SAT02.
2. The timeline is intentionally a two-level plot:
   - `coast` means \(u_T=0\), waiting at the 350 km parking shell for \(J_2\) RAAN drift.
   - `thrust` means the reduced Lyapunov raise is active.
3. The safety check is now an **SMA safety check**. The 200 km altitude limit is
   plotted as \(a_{\mathrm{safe}}=R_E+200\ \mathrm{km}\) on each transfer subplot,
   so overlapping transfers do not hide SAT01 or SAT02.

The coast time used for each satellite-slot pair is computed from

\[
\Delta\Omega_0+
(\dot\Omega_p-\dot\Omega_f)t_{\mathrm{coast}}
+\Delta\Omega_{\mathrm{raise}}\approx 0,
\]

so

\[
t_{\mathrm{coast}}=
-\frac{\Delta\Omega_0+\Delta\Omega_{\mathrm{raise}}}
{\dot\Omega_p-\dot\Omega_f}.
\]
