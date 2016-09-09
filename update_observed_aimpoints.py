#!/usr/bin/env python

"""
Trending application to update the ``OBSERVED_AIMPOINTS_FILE`` table (ascii ECSV format)
in place to reflect information about observed aimpoints, and in particular the delta
offset from the planned value.  This also makes a plot (default 6-months of data) for
inspection.
"""
import os
import glob
import functools
import shelve
import argparse

import numpy as np
from astropy.table import Table, vstack
from mica.starcheck import get_mp_dir
import Ska.Shell
import Ska.arc5gl
from kadi import events
from Chandra.Time import DateTime
import pyyaks.logger

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpld3
from Ska.Matplotlib import plot_cxctime


MP_ROOT = '/data/mpcrit1/mplogs'
OBSERVED_AIMPOINTS_FILE = 'observed_aimpoints.dat'

ciaoenv = Ska.Shell.getenv('source /soft/ciao/bin/ciao.sh')
ciaorun = functools.partial(Ska.Shell.bash, env=ciaoenv)
dmcoords_cmd = ['dmcoords', 'none',
                'asolfile=none',
                'detector="{detector}"',
                'fpsys="{fpsys}"',
                'opt=cel',
                'ra={ra_targ}',
                'dec={dec_targ}',
                'celfmt=deg',
                'ra_nom={ra_nom}',
                'dec_nom={dec_nom}',
                'roll_nom={roll_nom}',
                'ra_asp=")ra_nom"',
                'dec_asp=")dec_nom"',
                'roll_asp=")roll_nom"',
                'sim="{sim_x} 0 {sim_z}"',
                'displace="0 {dy} {dz} 0 0 0"',
                'verbose=0']
dmcoords_cmd = ' '.join(dmcoords_cmd)

opt = None

# Set up logging
loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(name='observed_aimpoints_mon', level=loglevel,
                                  format="%(asctime)s %(message)s")


class TooManyFilesError(ValueError):
    pass


class NoObsidError(ValueError):
    pass


def get_opt(args=None):
    parser = argparse.ArgumentParser(description='Plot aimpoint drift data '
                                     'from aspect solution files')
    parser.add_argument("--start",
                        help="Processing start date (default=NOW - 14 days)")
    parser.add_argument("--stop",
                        help="Processing stop date (default=NOW)")
    parser.add_argument("--data-root",
                        default=".",
                        help="Root directory for data files (default='.')")
    parser.add_argument("--lookback",
                        default=180,
                        help="Lookback time for plotting (days, default=180)")
    return parser.parse_args(args)


def get_evt_meta(obsid, detector):
    """
    Get event file metadata (FITS keywords) for ``obsid`` and ``detector`` and cache for later use.

    Returns a dict of key=value pairs, or None if there is no data in archive.
    """
    evts = shelve.open('event_meta.shelf')
    sobsid = str(obsid)
    if sobsid not in evts:
        det = 'hrc' if detector.startswith('HRC') else 'acis'
        arc5gl = Ska.arc5gl.Arc5gl()
        arc5gl.sendline('obsid={}'.format(obsid))
        arc5gl.sendline('get {}2'.format(det) + '{evt2}')
        del arc5gl

        files = glob.glob('{}f{}*_evt2.fits.gz'.format(det, obsid))
        if len(files) == 0:
            raise NoObsidError('No event file found for obsid {}'.format(obsid))
        if len(files) > 1:
            raise TooManyFilesError('Wrong number of files {}'.format(files))

        evt2 = Table.read(files[0], hdu=1)
        os.unlink(files[0])

        evts[sobsid] = {k.lower(): v for k, v in evt2.meta.items()}

    out = evts[sobsid]
    evts.close()

    return out


def dmcoords_chipx_chipy(keys, verbose=False):
    """
    Get the dmcoords-computed chipx and chipy for given event file
    header keyword params.  NOTE: the ``dy`` and ``dz`` inputs
    to dmcoords are flipped in sign from the ASOL values.  Generally the
    ASOL DY/DZ are positive and dmcoord input values are negative.  This
    sign flip is handled *here*, so input to this is ASOL DY/DZ.

    :param keys: dict of event file keywords
    """
    # See the absolute_pointing_uncertainty notebook in this repo for the
    # detailed derivation of this -15.5, 6.0 arcsec offset factor.  See the
    # cell below for the summary version.
    ciaorun('punlearn dmcoords')
    fpsys_map = {'HRC-I': 'HI1',
                 'HRC-S': 'HS2',
                 'ACIS': 'ACIS'}
    keys = {key.lower(): val for key, val in keys.items()}
    det = keys['detnam']
    keys['detector'] = (det if det.startswith('HRC') else 'ACIS')
    keys['dy'] = -keys['dy_avg']
    keys['dz'] = -keys['dz_avg']
    keys['fpsys'] = fpsys_map[keys['detector']]

    cmd = dmcoords_cmd.format(**keys)
    ciaorun(cmd)

    if verbose:
        print(cmd)
    return [float(x) for x in ciaorun('pget dmcoords chipx chipy chip_id')]


def get_mp_aimpoint_entry(obsid):
    """
    Return an astropy table Row for the corresponding row in the *_dynamical_offsets.txt
    table in the SOT MP backstop products for the given ``obsid``.

    :param obsid: obsid (int)
    :returns: astropy.table.Row
    """
    mp_dir, status, run_date = get_mp_dir(obsid)
    if status != 'ran':
        raise NoObsidError('obsid {} has not run yet, status={!r}'.format(obsid, status))
    mp_dir = MP_ROOT + mp_dir  # mp_dir looks like '/2016/AUG2916/oflsa/'

    files = glob.glob(os.path.join(mp_dir, 'output', '*_dynamical_offsets.txt'))
    if len(files) != 1:
        raise ValueError('found {} dynamical offsets files in {}'.format(len(files), mp_dir))

    aimpoints = Table.read(files[0], format='ascii.basic', guess=False)
    idxs = np.flatnonzero(aimpoints['obsid'] == obsid)
    if len(idxs) != 1:
        raise ValueError('found {} matches for obsid {} in {}'
                         .format(len(idxs), obsid, files[0]))
    return aimpoints[idxs[0]]


def get_observed_aimpoint_offset(obsid):
    """
    Compare the predicted CHIPX/Y values with planned using observed event file
    data on actual ACA alignment.
    """
    plan = get_mp_aimpoint_entry(obsid)
    detector = plan['detector']

    evt = get_evt_meta(obsid, detector)

    obs_chipx, obs_chipy, obs_chip_id = dmcoords_chipx_chipy(evt)

    arcsec_per_pixel = 0.13175 if detector.startswith('HRC') else 0.492
    arcsec_per_mm = 20.49  # arcsec/radian / focal length (mm)

    plan_chipx = plan['chipx']
    plan_chipy = plan['chipy']
    offset_y = plan['target_offset_y']
    offset_z = plan['target_offset_z']

    # Coordinates:
    # In cycle 18 POG figures 4.26 - 29.
    #
    # SIM_Z:
    # -offset_z + arcsec_per_mm * (sim_z_nom - evt['sim_z'])
    #
    # ACIS-S : -offset_y <=> SIM +Y <=> +chip_x, -offset_z <=> SIM +Z <=> +chip_y
    # ACIS-I : -offset_y <=> SIM +Y <=> -chip_y, -offset_z <=> SIM +Z <=> +chip_x
    # HRC-S : -offset_y <=> SIM +Y <=> -chip_y, -offset_z <=> SIM +Z <=> +chip_x
    # HRC-I : -offset_y <=> SIM +Y <=> +chip_x / sqrt(2) - chip_y / sqrt(2)
    #       : -offset_z <=> SIM +Z <=> -chip_x / sqrt(2) - chip_y / sqrt(2)
    #          (-offset_y - offset_z) / sqrt(2) <=> chip_y
    #          (-offset_y + offset_z) / sqrt(2) <=> chip_x

    sim_z_nom = {'ACIS-S': -190.1401,
                 'ACIS-I': -233.5874,
                 'HRC-S': 250.4660,
                 'HRC-I': 126.98298}

    sim_z_off = arcsec_per_mm * (sim_z_nom[detector] - evt['sim_z'])  # arcsec

    if detector == 'ACIS-S':
        delta_x = -offset_y
        delta_y = -offset_z + sim_z_off

    elif detector == 'ACIS-I':
        delta_x = -offset_z + sim_z_off
        delta_y = +offset_y

    elif detector == 'HRC-S':
        delta_x = -offset_z + sim_z_off
        delta_y = +offset_y

    elif detector == 'HRC-I':
        # What are the odds this is right?
        delta_x = (-offset_y + offset_z + sim_z_off) / np.sqrt(2)
        delta_y = (-offset_y - offset_z - sim_z_off) / np.sqrt(2)

    else:
        raise ValueError('illegal detector {}'.format(detector))

    plan_chipx += delta_x / arcsec_per_pixel
    plan_chipy += delta_y / arcsec_per_pixel

    dx = (plan_chipx - obs_chipx) * arcsec_per_pixel
    dy = (plan_chipy - obs_chipy) * arcsec_per_pixel

    out = {'dx': dx,
           'dy': dy,
           'dr': np.hypot(dx, dy),
           'obs_chipx': obs_chipx,
           'obs_chipy': obs_chipy
           }

    for colname in plan.colnames:
        out[colname] = plan[colname]

    return out


def update_observed_aimpoints():
    """
    Update the ``OBSERVED_AIMPOINTS_FILE`` table (ascii ECSV format) in
    place to reflect information about observed aimpoints, and in particular
    the delta offset from the planned value.
    """
    # Default is between NOW and NOW - 14 days
    start = DateTime(opt.start) - (14 if opt.start is None else 0)
    stop = DateTime(opt.stop)

    # Get science obsids
    obsids = [evt.obsid for evt in events.obsids.filter(start, stop)
              if evt.obsid < 40000]

    # Read in existing file if it exists and make a set of already-processed obsids
    filename = os.path.join(opt.data_root, OBSERVED_AIMPOINTS_FILE)
    if os.path.exists(filename):
        logger.info('Reading {}'.format(filename))
        dat_old = Table.read(filename, format='ascii.ecsv', guess=False)
        processed_obsids = set(dat_old['obsid'])
    else:
        dat_old = None
        processed_obsids = set()

    rows = []
    for obsid in obsids:
        if obsid in processed_obsids:
            logger.info('Skipping obsid {}: already processed'.format(obsid))
            continue

        try:
            vals = get_observed_aimpoint_offset(obsid)
        except NoObsidError:  # not yet in archive
            logger.info('Skipping obsid {}: not in archive yet'.format(obsid))
            continue
        except Exception as err:
            logger.info('ERROR: {}'.format(err))
            continue

        logger.info('Obsid={obsid:5d} detector={detector:6s} '
                    'chipx={chipx:.1f} chipy={chipy:.1f} dr={dr:.1f}'
                    .format(**vals))
        rows.append(vals)

    if rows:
        dat = Table(rows=rows, names=sorted(rows[0]))
        if dat_old is not None:
            dat = vstack([dat_old, dat])
        logger.info('Writing {}'.format(filename))

        for name in 'dr dx dy obs_chipx obs_chipy'.split():
            dat[name].format = '.2f'

        dat.sort('mean_date')
        dat.write(filename, format='ascii.ecsv')
    else:
        dat = dat_old

    return dat


def plot_observed_aimpoints(obs_aimpoints):
    """
    Make png and html (mpld3) plot of data in the ``obs_aimpoints`` table.
    """
    plt.close(1)
    fig = plt.figure(1, figsize=(8, 4))

    dates = DateTime(obs_aimpoints['mean_date'])
    years = dates.frac_year
    times = dates.secs
    ok = years > np.max(years) - opt.lookback / 365.25
    obs_aimpoints = obs_aimpoints[ok]
    times = times[ok]

    plot_cxctime(times, obs_aimpoints['dx'], 'ob', label='CHIPX')
    plot_cxctime(times, obs_aimpoints['dy'], 'or', label='CHIPY')
    plt.grid()
    ymax = max(5, np.max(np.abs(obs_aimpoints['dx'])), np.max(np.abs(obs_aimpoints['dy'])))
    plt.ylim(-ymax, ymax)
    plt.ylabel('Offset (arcsec)')
    plt.title('Observed aimpoint offsets')

    plt.legend(loc='upper left', fontsize='small', title='')

    outroot = os.path.join(opt.data_root, 'observed_aimpoints')
    logger.info('Writing plot files {}.png,html'.format(outroot))
    mpld3.plugins.connect(fig, mpld3.plugins.MousePosition(fmt='.1f'))
    mpld3.save_html(fig, outroot + '.html')
    fig.patch.set_visible(False)
    plt.savefig(outroot + '.png', frameon=False)


def main():
    global opt
    opt = get_opt()

    obs_aimpoints = update_observed_aimpoints()
    plot_observed_aimpoints(obs_aimpoints)


if __name__ == '__main__':
    main()
