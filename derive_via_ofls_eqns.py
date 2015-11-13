# -*- coding: utf-8 -*-

from __future__ import print_function
from Quaternion import Quat

# Confirm that Quat.transform gives the rotation matrix that
# is specified in the OFLS attitudes and offsets definition document.

q_ra = Quat([20, 0, 0])
print('q_ra')
print(q_ra.transform)
# Confirm that this is R3(−α) = R3(-20 deg): YES

q_dec = Quat([0, 40, 0])
print('q_dec')
print(q_dec.transform)
# Confirm that this is R2(δ) = R2(40 deg): YES

q_roll = Quat([0, 0, 60])
print('q_roll')
print(q_roll.transform)
# Confirm that this is R1(-phi) = R1(-60 deg): YES

q_ra_dec_roll = Quat([20, 40, 60])

# From OFLS eqns:
#
#   V_I = R3(−α) R2(δ) R1(-φ) X = q.transform . X
#
# This rotates the X-axis vector in ECI.
#
M = q_ra.transform .dot (q_dec.transform) .dot (q_roll.transform)
print('M\n', M)
print('q_ra_dec_roll.transform\n', q_ra_dec_roll.transform)
print((q_ra * q_dec * q_roll).transform)

print(np.all(M == q_ra_dec_roll.transform))

# SO...
#
# q([RA, Dec, Roll]).transform is the rotation matrix M
# that rotates the unit X vector [1, 0, 0] to be coincident
# with the unit vector pointing at RA, Dec.

# Equivalently it rotates the ECI frame to be coincident with the
# spacecraft frame at the attitude defined by q.
#
# This corresponds to active rotation of a vector in ECI.
#
# The inverse of q.transform is the matrix which gives the
# components of an ECI vector in the spacecraft frame
# representation specified by q.   In other words:
#   V_q = M_q . V_I = q.transform.inv() . V_I


# The inertial to q-frame matrix is for q = Quat([alpha, delta, theta]):
#
# q.transform.inv() = M_q(φT)= R1(φT) R2(−δT) R3(αT)

# Inertial to ACA transform matrix = q_aca.transform.inv() is:
#
# M_AI(φT) = M_AS R3(−ΔY) R2(ΔZ) R1(φT) R2(−δT) R3(αT)
#          = M_AS R3(−ΔY) R2(ΔZ) q_targ.transform.inv()
# q_aca.transform.inv() = M_AS R3(−ΔY) R2(ΔZ) q_targ.transform.inv()
# q_aca.transform.inv() = M_AS . Quat([ΔY,ΔZ]).transform .  q_targ.transform.inv()

# Finally:
# q_aca.transform = q_targ.transform . q_off.transform.inv() . M_AS.inv()

# ------------

# Original equations in calc_si_align specified:
#  q_si_align = Quat(si_align)
#  q_targ = Quat([ra_targ, dec_targ, roll_targ])
#  q_off = Quat([-y_off, -z_off, 0])
#  q_pcad = q_targ * q_off * q_si_align

# si_align was transposed, which is the same as inverse for a rotation matrix.

# q_off was defined with -y_off and -z_off, which is the same (to first order)
# as the inverse of Quat([y_off, z_off, 0])

# Therefore we can write the equations *correctly* as:

# With ODB_SI_ALIGN correctly defined in C order, which effectively
# inverts the original q_si_align:
#  q_si_align = Quat(ODB_SI_ALIGN)
#  q_targ = Quat([ra_targ, dec_targ, roll_targ])
#  q_off = Quat([y_off, z_off, 0])
#  q_pcad = q_targ * q_off.inv() * q_si_align.inv()
