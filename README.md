# SINQCharta

A collection of useful functions for binning, adding and subtracting
neutron scattering scans from SINQ instruments.

Running the following will display this `README` as well as the list of
exported functions of the module.

```julia
julia> ?
help?> SINQCharta
```

Considered data formats include:

  - ILL-formatted files from the spectrometers `TASP` and `EIGER`.
  - General ASCII files from `ZEBRA` in point-detector mode.
  - ASCII files from `MPMS\PPMS` measurement devices.