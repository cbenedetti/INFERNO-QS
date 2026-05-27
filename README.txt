*** Copyright Notice ***

INtegrated Fluid & paRticles simulatioN cOde (INF&RNO) Copyright (c) 2026,
The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

================================================================================
                  LASER-DRIVEN PLASMA WAKEFIELD ACCELERATION
================================================================================

OVERVIEW
--------
This code simulates the propagation of a laser pulse through a plasma channel
using the quasi-static module of the INF&RNO particle-in-cell (PIC) code.
The wakefield equations are solved self-consistently with the laser envelope evolution.
An externally injected Gaussian electron bunch is tracked through the plasma wake.

The simulation supports user-defined plasma density profiles with entrance and
exit ramps, a parabolic transverse channel for laser guiding, and full
diagnostics output at user-specified intervals.


--------------------------------------------------------------------------------
REPOSITORY CONTENTS
--------------------------------------------------------------------------------

  main.c                                  - Main simulation file
  source_QS_laserIndependentGrid.c - Quasi-static solver and field routines
  source_beamInitialDistribution.c - Beam generation and particle push routines
  README.txt                      - This file


--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------

  - GCC or compatible C compiler (C99 standard)
  - Standard C math library (libm)
  - source_QS_laserIndependentGrid.c must be present in the same directory
  - source_beamInitialDistribution.c must be present in the same directory

--------------------------------------------------------------------------------
COMPILATION
--------------------------------------------------------------------------------

  gcc main.c -O3 -std=c99  -lm

--------------------------------------------------------------------------------
EXECUTION
--------------------------------------------------------------------------------

 ./a.out

All the output files (see below) are produced in the local directory ./

--------------------------------------------------------------------------------
UNITS
--------------------------------------------------------------------------------

  All quantities use normalized plasma units unless otherwise noted.

    Length      normalized to the plasma skin depth  kp^-1 = c / omega_p
    Fields      normalized to  E0 = m_e c omega_p / e
    Vector pot. normalized to  a = e A / m_e c^2
    Energy      expressed as the Lorentz factor  gamma = E / m_e c^2

--------------------------------------------------------------------------------
PHYSICS
--------------------------------------------------------------------------------

  Plasma Density Profile
  ----------------------
  The plasma is defined along the propagation axis z as follows:

    [ vacuum ][ ramp in ][ flat top ][ ramp out ][ vacuum ]

  The transverse profile follows a parabolic channel:

    n(z, r) = [1 + 4r^2 / matched_radius^4] x longitudinal_envelope(z)

  This channel is designed to guide the laser by matched propagation.

 The function defining the plasma must be defined by the user (see the function "plasma" in  main.c)
 
  Laser Pulse
  -----------
  The laser is modeled with an initially bi-Gaussian envelope and evolves self-consistently in the plasma.
  The key dimensionless parameters are the normalized vector potential a0, pulse length kpL, spot size kpW, and the
  wavenumber ratio k0_kp.

  Injected Beam
  -------------
  A Gaussian electron bunch is injected externally. The distribution is generated with the routine generateBeam_GaussianXGaussian,
  which initializes the bunch with prescribed position, size, emittance, and energy.

--------------------------------------------------------------------------------
KEY PARAMETERS (the numerical values provided are the ones in the current main.c_
--------------------------------------------------------------------------------

  Quasi-Static Solver
  -------------------
  iter          300       Maximum iterations per propagation step
  err           1e-4      Convergence tolerance
  alpha_A       1.0       Relaxation parameter 
  alpha_B       0.98      Relaxation parameter 
  psi_min_0    -0.99      Initial minimum pseudo-potential
  d_psi_min     0.01      Pseudo-potential increment

  Plasma
  ------
  z_plasma_start    6.0      Start of plasma  [kp^-1]
  L_ramp_entrance   100.0    Entrance ramp length  [kp^-1]
  L_plasma          1000.0   Flat-top plasma length  [kp^-1]
  L_ramp_exit       100.0    Exit ramp length  [kp^-1]
  matched_radius    3.0      Matched channel radius  [kp^-1]

  Laser (bi-Gaussian laser, specified by the function "gaussian_envelope")
  -----
  a0            2.0       Normalized peak vector potential
  kpL           1.5       Normalized pulse length
  kpW           3.0       Normalized spot size
  k0_kp         100.0     Ratio of laser to plasma wavenumber
  kpz0          0.0       Initial longitudinal focus position
  kpzf          0.0       Final longitudinal focus position

  Injected Bunch
  --------------
  Npb           20000     Number of macro-particles
  kpz0b        -5.0       Initial longitudinal centroid position  [kp^-1]
  kpsigmazb     0.3       Longitudinal RMS size  [kp^-1]
  kpsigmaxb     0.3       Transverse RMS size  [kp^-1]
  kpepsib       0.5       Normalized emittance  [kp^-1]
  gammab        2000      Initial Lorentz factor
  dgammab       0.0       Initial energy spread
  nb_over_n0b   1.0       Bunch-to-plasma density ratio
  selfConsb=1 indicates the bunch is self-consistent. If selfConsb=NO the bunch is passive (i.e., it does not generate a wake).

  Computational Grid
  ------------------
  zmin         -8.0       Rear of co-moving window  [kp^-1]
  zmax          6.0       Front of co-moving window  [kp^-1]
  dz            1/30      Longitudinal resolution, wakefield grid  [kp^-1]
  dzl           1/200     Longitudinal resolution, laser grid  [kp^-1]
  rmax          16.0      Transverse box size  [kp^-1]
  dr            1/20      Transverse resolution, wakefield grid  [kp^-1]
  drl           1/10      Transverse resolution, laser grid  [kp^-1]
  ds            Zr/300    Propagation step size  [kp^-1] (expressed as a fraction of the laser's Rayleigh length)
  Nppc          2         Particles per cell

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------

  All output is written to DIRECTORY_OUTPUT (default: current directory "./" ).

  densityProfile    On-axis plasma density n(z,r=0)/n0 vs. propagation distance [ASCII]
  a0_max             Peak normalized vector potential vs. propagation distance [ASCII]
  laserEnergy       Laser energy (noramalized to the initial value) vs. propagation distance [ASCII]
  fieldData           Dump of all field data at snapshot steps [BINARY] 
  cutEz                 Lineout of the longitudinal electric field Ez on axis at snapshot steps [ASCII]
  bunch               Dump of the full beam phase space at snapshot steps [ASCII]

  Snapshots are saved at 131 evenly spaced intervals of 10 kp^-1, starting
  from s = 0. Distrbution of the oputputs can be changed by the user by modifying the array dumpList[]

--------------------------------------------------------------------------------
SIMULATION LOOP
--------------------------------------------------------------------------------

  Each propagation step s performs the following sequence:

    1. Deposit beam charge onto the grid (linear shape function)
    2. Solve the quasi-static wake equations
    3. Record peak laser amplitude and total laser energy
    4. Write field and particle snapshots if the step matches the dump list
    5. Evolve the laser pulse envelope
    6. Push beam particles through the wake fields


================================================================================
```
