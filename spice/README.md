### SPICE kernel folder

The SPICE kernels should be placed in this folder for use in the program subroutines. These should be loadable
through a metakernel `metakernel.tm`, which should contain the locations and names of the kernels.

The kernels required for program run-time are:

* A leapseconds file (`naif0012.tls` recommended)
* A planetary ephemeris file (`de430.bsp` recommended)
* A planetary masses file (`de-403-masses.tpc` recommended)
* A planetary PCK file (`pck00010.tpc` recommended)
* Any additional kernels for non-standard solar system bodies.

The metakernel structure should look like this:

```
   KPL/MK
   \begindata
 
   KERNELS_TO_LOAD = ( 'leapseconds_file.tls',
                       'planetary_ephemeris_file.bsp',
                       'planetary_masses_file.tpc',
                       'planetary_pck_file.tpc')
 
   \begintext
```