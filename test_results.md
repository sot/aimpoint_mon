# Testing of automated ODB_SI_ALIGN updates

### Test sign of aimpoint offset

The [aimpoint trending page](http://cxc.cfa.harvard.edu/mta/ASPECT/aimpoint_mon/) shows
the daily update to the aimpoint offsets for ACIS-S and ACIS-I.  This is computed
within the `plot_aimpoint.py` module.

The aimpoint offsets specify values which will move the mean observed aimpoint so
that it matches the POG aimpoints for the ACIS SIs.  (HRC is not adjusted because
there are no pointing sensitive observations with HRC).

The signs of the values are crucial, but performing an end-to-end test is not practical.
The following three independent analyses have been done:

- T. Aldcroft: used the aimpoint drift plots in Chapter 4 of the POG, which include arrows
  that indicate the direction of an applied DY or DZ offset.
- A. Tennant: used ObsVis to test with Obsid 16361 (April Crab observation with a
  large offset).
- J. McDowell: used DS tool `dmcoords` to independently confirm the derived relationship
  between *aspect solution* DY, DY and CHIPX, CHIPY.
- J. Connelly used ObsVis in a different way from Tennant to confirm the sign and magnitude of
offsets

**Test pass: Yes**

### Test computation of updated ODB_SI_ALIGN values

**Strategy**: use a real observation with non-trivial and non-identical
OR-specified target offsets.  This provides:

- Target RA, Dec (science target, not OBC target)
- Target offset DY, DZ
- PCAD pointing attitude that gets uplinked via command loads.
  The PCAD attitude is at the target attitude with the ODB_SI_ALIGN
  transform applied and the specified target offset applied.

The aimpoint alignment offset from http://cxc.cfa.harvard.edu/mta/ASPECT/aimpoint_mon/
is defined as an ADDITIVE offset to the Target offset.  Specifically, define:

```
DY_total = DY_target + dy_align
DZ_total = DZ_target + dz_align
```

These `DY_total`, `DZ_total` offsets *could* be used in place of the OR-specified values
in the planning process to correct for aimpoint drift.  But in fact the aimpoint drift
is being corrected via the ODB_SI_ALIGN matrix.

A function `calc_si_align(dy_align, dz_align)` computes the new SI_ALIGN matrix for
given dy_align, dz_align values.

To test that this is giving correct values (and the signs in particular), *pretend*
that the required aimpoint alignment offset is exactly the Target offset from
obsid 13928.  In this case we compute a new SI_ALIGN with that offset and then zero out
the target offset.  By the additive property above the final PCAD attitude must be the
same as the originally scheduled (and uplinked) attitude.

The relevant OR and pointing attitude for obsid 13928 are:

```
Observation request:
 ID=13928,TARGET=(191.321250,27.125556,{Haro 9}),DURATION=(17000.000000),
 PRIORITY=9,SI=ACIS-S,GRATING=NONE,SI_MODE=TE_0045A,ACA_MODE=DEFAULT,
 TARGET_OFFSET=(0.002500,-0.004167),
 DITHER=(ON,0.002222,0.360000,0.000000,0.002222,0.509100,0.000000),
 SEGMENT=(1,15300.000000),PRECEDING=(13632),MIN_ACQ=1,MIN_GUIDE=1

PCAD:  As-planned pointing from starcheck
  Q1,Q2,Q3,Q4: -0.18142595  -0.37811633  -0.89077416  0.17502588
```

This test is encoded in the function `test_si_align_with_data()` within the
`calc_si_align` module.

**Test pass: Yes**

### Test updated Characteristics file in OFLS

A new characteristics file `CHARACTERIS_09OCT15` was created on 2015-Oct-9 using the
function call:

```
update_characteristics('CHARACTERIS_12MAR15',
                        dy_acis_i=-5, dz_acis_i=10,
                        dy_acis_s=15, dz_acis_s=-20)
```

The OFLS was used to generate products using both the new `CHARACTERIS_09OCT15` and the
current baseline `CHARACTERIS_12MAR15` characteristics files.  The output products are
found in the `test_ofls` directory (available on request but not part of the git
repository).  The script `test_ofls/compare_attitudes.py` compares each of the final
attitudes within the maneuver summary file.

The listing below shows the script output with additional annotation indicating the
science instrument.  In all cases the expected values of `DY` and `DZ` are seen.  The roll
offsets are as large as 41 arcsec.  The source of this is not precisely understood, but
likely stems from using the original PCAD commanded roll as the target roll.  This offset
is comparable to typical variations in roll during normal point mode, and corresponds to a
0.1 arcsec displacement at 8 arcmin off-axis (i.e. it is negligible).

```
Obsid 17503:  DY= -5.00  DZ= 10.00  droll=  2.20 arcsec   ACIS-I
Obsid P3201:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3202:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3203:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3204:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid 17244:  DY= 15.00  DZ=-20.00  droll=  0.01 arcsec   ACIS-S
Obsid 17283:  DY= -5.00  DZ= 10.00  droll=  4.26 arcsec   ACIS-I
Obsid 17284:  DY= -5.00  DZ= 10.00  droll=  4.22 arcsec   ACIS-I
Obsid 17285:  DY= -5.00  DZ= 10.00  droll=  4.19 arcsec   ACIS-I
Obsid 17286:  DY= -5.00  DZ= 10.00  droll=  4.19 arcsec   ACIS-I
Obsid 17008:  DY= -5.00  DZ= 10.00  droll=  6.83 arcsec   ACIS-I
Obsid 17027:  DY= 15.00  DZ=-19.99  droll= 13.40 arcsec   ACIS-S
Obsid 16648:  DY= -5.00  DZ= 10.00  droll= -2.99 arcsec   ACIS-I
Obsid P3301:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3302:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid P3303:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid 17011:  DY= -5.00  DZ= 10.00  droll=  6.93 arcsec   ACIS-I
Obsid 17014:  DY= -5.00  DZ= 10.00  droll=  0.00 arcsec   ACIS-I
Obsid 17119:  DY= 15.00  DZ=-20.00  droll= -6.31 arcsec   ACIS-S
Obsid 16670:  DY= -5.00  DZ= 10.00  droll=  2.55 arcsec   ACIS-I
Obsid 17552:  DY= -5.00  DZ= 10.00  droll=  4.15 arcsec   ACIS-I
Obsid P3401:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3402:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid P3403:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3404:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid 17527:  DY= -5.00  DZ= 10.00  droll= -0.02 arcsec   ACIS-I
Obsid 17012:  DY= -5.00  DZ= 10.00  droll= -0.00 arcsec   ACIS-I
Obsid 18188:  DY= 15.00  DZ=-20.00  droll=  0.02 arcsec   ACIS-S
Obsid 18669:  DY= 15.00  DZ=-20.00  droll=  0.01 arcsec   ACIS-S
Obsid 16763:  DY= 15.00  DZ=-20.00  droll=  7.89 arcsec   ACIS-S
Obsid P3501:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3502:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid P3503:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid 18682:  DY= -5.00  DZ= 10.00  droll= -0.02 arcsec   ACIS-I
Obsid 17446:  DY= -5.00  DZ= 10.00  droll=  7.42 arcsec   ACIS-I
Obsid 18684:  DY= 15.00  DZ=-20.00  droll=  0.02 arcsec   ACIS-S
Obsid 16825:  DY= -5.00  DZ= 10.00  droll= -1.86 arcsec   ACIS-I
Obsid 18681:  DY= 15.00  DZ=-20.00  droll=  0.01 arcsec   ACIS-S
Obsid 18641:  DY= -5.00  DZ= 10.00  droll=  1.15 arcsec   ACIS-I
Obsid 17087:  DY= 15.01  DZ=-20.02  droll=-41.04 arcsec   ACIS-S
Obsid P3601:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3602:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid P3603:  DY= -0.00  DZ=  0.00  droll=  0.00 arcsec     --
Obsid P3604:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
Obsid 17529:  DY= -5.00  DZ= 10.00  droll= -0.01 arcsec   ACIS-I
Obsid 17013:  DY= -5.00  DZ= 10.00  droll=  7.41 arcsec   ACIS-I
Obsid 16662:  DY= 15.00  DZ=-19.99  droll= 13.59 arcsec   ACIS-S
Obsid 16827:  DY= -5.00  DZ= 10.00  droll=  0.00 arcsec   ACIS-I
Obsid 18685:  DY= 15.00  DZ=-20.00  droll=  0.02 arcsec   ACIS-S
Obsid 16945:  DY= 15.00  DZ=-20.00  droll=  0.00 arcsec   ACIS-S
Obsid 17097:  DY= 15.00  DZ=-20.00  droll=  0.01 arcsec   ACIS-S
Obsid 18683:  DY= -5.00  DZ= 10.00  droll= -0.01 arcsec   ACIS-I
Obsid P3701:  DY= -0.00  DZ= -0.00  droll=  0.00 arcsec     --
```

**Test pass: Yes**

### Regression test after updates to calc_si_align and plot_aimpoints

Following comments from community, calc_si_align code was updated
so that the ODB_SI_ALIGN matrix was transposed.  This makes the
matrix correspond to the OFLS matrix.

Using plot_aimpoints.py with the AIMPOINT_JUMPS variable stubbed
out, a new local info.json was created.  This was then used
to compute the CHARACTERIS_03NOV15 files.

Confimed that diffs from CHARACTERIS_12OCT15 were as expected
(basically just negligible numerical diffs).
