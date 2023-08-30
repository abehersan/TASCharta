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

  - ILL-formatted files from `TASP` and `EIGER`.
  - General ASCII files from `ZEBRA` in point-detector mode.
  - HDF files from `DMC` using functions afforded by `DMCpy`.
  - HDF files from `CAMEA` using functions afforded by `MJOLNIR`. 
