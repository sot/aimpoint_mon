import numpy as np
from Quaternion import Quat

# Baseline ODB_SI_ALIGN used prior to Nov. 2015.  This is expressed to
# 5 digits of precision and is not an orthonormal rotation matrix.

ODB_SI_ALIGN = np.array([[1.0, 3.3742E-4, 2.7344E-4],
                         [-3.3742E-4, 1.0, 0.0],
                         [-2.7344E-4, 0.0, 1.0]])

# Fix the transform matrix by converting to a normalized Quaternion
# and then back.  This gives:
# array([[  9.99999906e-01,   3.37420000e-04,   2.73440000e-04],
#        [ -3.37420000e-04,   9.99999943e-01,  -4.61320624e-08],
#        [ -2.73440000e-04,  -4.61320624e-08,   9.99999963e-01]])

ODB_SI_ALIGN = Quat(Quat(ODB_SI_ALIGN).q).transform

ARCSEC2RAD = np.pi / (3600 * 180)


def test_si_align_with_data():
    """
    Test an actual observation planned with a target offset in Y and Z:

    Strategy: use a real observation with non-trivial and non-identical
    OR-specified target offsets.  This provides:

      - Target RA, Dec (science target, not OBC target)
      - Target offset DY, DZ
      - PCAD pointing attitude that gets uplinked via command loads.
        The PCAD attitude is at the target attitude with the ODB_SI_ALIGN
        transform applied and the specified target offset applied.

    The aimpoint alignment offset from http://cxc.cfa.harvard.edu/mta/ASPECT/aimpoint_mon/
    is defined as an ADDITIVE offset to the Target offset.  Specifically, define:

     DY_total = DY_target + dy_align
     DZ_total = DZ_target + dz_align

    These DY_total, DZ_total offsets *could* be used in place of the OR-specified values
    in the planning process to correct for aimpoint drift.  But in fact the aimpoint drift
    is being corrected via the ODB_SI_ALIGN matrix.

    A function calc_si_align(dy_align, dz_align) computes the new SI_ALIGN matrix for
    given dy_align, dz_align values.

    To test that this is giving correct values (and the signs in particular), *pretend*
    that the required aimpoint alignment offset is exactly the Target offset from obsid
    13928.  In this case we compute a new SI_ALIGN with that offset and then zero out
    the target offset.  By the additive property above the final PCAD attitude must be the
    same as the originally scheduled (and uplinked) attitude.

    ----------------------

    Observation request:
     ID=13928,TARGET=(191.321250,27.125556,{Haro 9}),DURATION=(17000.000000),
     PRIORITY=9,SI=ACIS-S,GRATING=NONE,SI_MODE=TE_0045A,ACA_MODE=DEFAULT,
     TARGET_OFFSET=(0.002500,-0.004167),
     DITHER=(ON,0.002222,0.360000,0.000000,0.002222,0.509100,0.000000),
     SEGMENT=(1,15300.000000),PRECEDING=(13632),MIN_ACQ=1,MIN_GUIDE=1

    PCAD:  As-planned pointing from starcheck
      Q1,Q2,Q3,Q4: -0.18142595  -0.37811633  -0.89077416  0.17502588
    """
    # Attitude quaternion for the as-run PCAD attitude
    q_pcad = Quat([-0.18142595, -0.37811633, -0.89077416, 0.17502588])

    # Target coordinates and quaternion, using the PCAD roll
    ra_targ, dec_targ = 191.321250, 27.125556
    q_targ = Quat([ra_targ, dec_targ, q_pcad.roll])

    # Offsets from OR (Target DY, DZ) in degrees
    dy_offset, dz_offset = 0.002500, -0.004167

    # Sanity check that offset between as-planned PCAD and target are as expected.  This
    # uses the baseline ODB_SI_ALIGN.  The calc_offsets() function has been independently
    # validated for planning purposes.
    dy_nom_align, dz_nom_align = calc_offsets(ra_targ, dec_targ,
                                              q_pcad.ra, q_pcad.dec, q_pcad.roll)
    assert np.allclose(dy_nom_align, dy_offset, atol=3e-6)  # 0.01 arcsec
    assert np.allclose(dz_nom_align, dz_offset, atol=3e-6)

    # Test calc_pcad_from_targ() for baseline SI_ALIGN
    q_pcad_calc = calc_pcad_from_targ(ra_targ, dec_targ, q_pcad.roll, dy_offset, dz_offset)
    dq = q_pcad.inv() * q_pcad_calc
    assert np.allclose(dq.q[0] * 2, 0, atol=36 * ARCSEC2RAD)
    assert np.allclose(dq.q[1] * 2, 0, atol=0.01 * ARCSEC2RAD)
    assert np.allclose(dq.q[2] * 2, 0, atol=0.01 * ARCSEC2RAD)

    # Compute PCAD attitude using a new SI ALIGN values based on OR-list TARGET_OFFSET
    # values and with no target offset applied.  This must match the as-planned
    # attitude using the baseline ODB_SI_ALIGN and WITH a target offset applied.
    si_align = calc_si_align(dy_offset, dz_offset, check_consistency=True)
    q_pcad_align = q_targ * Quat(si_align)

    dy_align, dz_align = calc_offsets(ra_targ, dec_targ, *q_pcad_align.equatorial)
    dq = q_pcad.inv() * q_pcad_align

    # Check delta quaternion (error) is less than 0.01 arcsec
    # in pitch and yaw, 36 arcsec in roll.
    assert np.allclose(dq.q[0] * 2, 0, atol=36 * ARCSEC2RAD)
    assert np.allclose(dq.q[1] * 2, 0, atol=0.01 * ARCSEC2RAD)
    assert np.allclose(dq.q[2] * 2, 0, atol=0.01 * ARCSEC2RAD)

    # Check offsets match to better than 0.01 arcsec
    assert np.allclose(dy_align, dy_offset, atol=3e-6)
    assert np.allclose(dz_align, dz_offset, atol=3e-6)


def calc_si_align(dy_align=0.0, dz_align=0.0, check_consistency=True):
    """
    Compute the SI alignment matrix that includes an alignment offset
    (dy_align, dz_align) from the baseline ODB_SI_ALIGN.

    :param dy_align: Aimpoint alignment offset in Dy (deg)
    :param dz_align: Aimpoint alignment offset in Dz (deg)
    :param check_consistency: double check the output by essentially inverting it

    :rtype: 3x3 alignment matrix
    """
    # Define new SI_ALIGN as a delta quaternion update to existing.
    # The minus signs here were empirically determined to pass tests.
    # This is the only place where signs are chosen without a clear justification.
    si_align = Quat([-dy_align, -dz_align, 0.0]) * Quat(ODB_SI_ALIGN)

    # Get the 3x3 transform matrix corresponding to the si_align quaternion
    out = si_align.transform

    # Double check that the si_align corresponds to the y and z offsets of the target in
    # that frame, in degrees.
    if check_consistency:
        q_pcad = Quat([1, 2, 3])  # Arbitrary attitude
        q_targ = q_pcad * si_align.inv()
        y_off, z_off = calc_offsets(q_targ.ra, q_targ.dec, q_pcad.ra, q_pcad.dec, q_pcad.roll)

        # Offset as expected?
        assert np.allclose(dy_align, y_off, atol=1e-6)
        assert np.allclose(dz_align, z_off, atol=1e-6)

        # Orthonormal?
        assert np.allclose(out.T.dot(out), np.eye(3))

    return out


def radec2yz(ra, dec, q):
    """
    Given target ``ra`` and ``dec`` and pointing quaternion ``q``, return Y and Z offset.
    The input ``ra`` and ``dec`` values can be 1-d arrays in which case the output ``dy``
    and ``dz`` will be corresponding arrays of the same length.

    :param ra: Target right Ascension (degrees)
    :param dec: Target declination (degrees)
    :param q: Pointing quaternion
    :rtype: list dy, dz (degrees)
    """
    r = np.radians(ra)
    d = np.radians(dec)
    eci = np.array([np.cos(r) * np.cos(d), np.sin(r) * np.cos(d), np.sin(d)])

    d_aca = np.dot(q.transform.transpose(), eci)

    dy = np.degrees(np.arctan2(d_aca[1], d_aca[0]))
    dz = np.degrees(np.arctan2(d_aca[2], d_aca[0]))

    return dy, dz


def calc_offsets(ra_targ, dec_targ, ra_pcad, dec_pcad, roll_pcad):
    """
    Calculates required Y and Z offsets (deg) required from a target to
    arrive at the desired PCAD pointing.

    :param ra_targ: RA of science target from OR/Obscat
    :param dec_targ: Dec of science target from OR/Obscat
    :param ra_pcad: RA of desired PCAD Pointing
    :param dec_pcad: Dec of desired PCAD Pointing
    :param roll_pcad: Roll of desired PCAD Pointing

    :rtype: tuple (y_off, z_off) arcmins
    """
    # Convert si_align transform matrix into a Quaternion
    si_align = Quat(ODB_SI_ALIGN)

    # Pointing quaternion
    q_pcad = Quat([ra_pcad, dec_pcad, roll_pcad])

    # Pointing quaternion of nominal HRMA frame after adjusting for the alignment offset.
    # The sense of si_align is that q_pcad = q_hrma * si_align, where si_align is
    # effectively a delta quaternion.
    q_hrma = q_pcad * si_align.inv()

    # the y and z offsets of the target in that frame, in degrees
    y_off, z_off = radec2yz(ra_targ, dec_targ, q_hrma)

    return y_off, z_off


def calc_pcad_from_targ(ra_targ, dec_targ, roll_targ, y_off, z_off, si_align=ODB_SI_ALIGN):
    """
    Calculates required Y and Z offsets (deg) required from a target to
    arrive at the desired PCAD pointing.

    :param ra_targ: RA of science target from OR/Obscat
    :param dec_targ: Dec of science target from OR/Obscat
    :param roll_targ: Roll of science target attitude
    :param y_off: Y offset (deg, sign per OR-list convention)
    :param z_off: Z offset (deg,sign per OR-list convention)
    :param si_align: SI ALIGN matrix (default=ODB_SI_ALIGN)

    :rtype: q_pcad (Quat)
    """
    q_si_align = Quat(si_align)
    q_targ = Quat([ra_targ, dec_targ, roll_targ])
    q_off = Quat([-y_off, -z_off, 0])
    q_pcad = q_targ * q_off * q_si_align

    return q_pcad
