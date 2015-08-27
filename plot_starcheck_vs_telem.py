#!/usr/bin/env python
"""
Plot positions of fid lights from starcheck versus the actual seen in telemetry.
::

  Usage: plot_starcheck_vs_telem.py [-h] [--start START] [--stop STOP]
                                    [--out OUT]

  Commanded vs. telemetry fid positions

  optional arguments:
    -h, --help     show this help message and exit
    --start START  Start date (default=NOW - 90 days)
    --stop STOP    Stop date (default=NOW)
    --out OUT      Output plot file
"""

import argparse
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astropy import table
from astropy.table import Table

from pyyaks.logger import get_logger
from kadi import events
import Ska.DBI
from Ska.Matplotlib import plot_cxctime
from Ska.engarchive import fetch_eng as fetch
from Chandra.Time import DateTime

logger = get_logger()


def get_opt(args):
    parser = argparse.ArgumentParser(description='Commanded vs. telemetry fid positions')
    parser.add_argument('--start',
                        type=str,
                        help='Start date (default=NOW - 90 days)')
    parser.add_argument('--stop',
                        type=str,
                        help='Stop date (default=NOW)')
    parser.add_argument('--out',
                        type=str,
                        help='Output plot file')

    opt = parser.parse_args(args)
    return opt


def get_fids_starcheck(dwells):
    """
    Get fid starcheck values corresponding to ``dwells`` (keyed by obsid)

    :returns: dict of tables keyed by obsid
    """

    logger.info('Getting fids from sybase starcheck table')
    fids_starcheck = {}
    with Ska.DBI.DBI(server='sybase', dbi='sybase', user='aca_read') as db:
        for obsid, dwell in dwells.items():
            # Only accept the first dwell of science observations
            if obsid in fids_starcheck:
                logger.info('Skipping obsid {} already in fids_starcheck'.format(obsid))
                continue

            fids = db.fetchall('SELECT obsid, slot, yang, zang FROM starcheck_catalog '
                               'WHERE obsid={} AND type="FID"'
                               .format(obsid))
            if len(fids) > 0:
                fids_starcheck[obsid] = Table(fids)
            else:
                logger.info('No fids found for obsid {} in starcheck_catalog'.format(obsid))

    return fids_starcheck


def get_dwells_with_fids(start, stop):
    """
    Get telemetry yag, zag values for each fid in ``fids_starcheck`` at the start
    of the corresponding ``dwells``.

    :returns: dict of tables keyed by obsid (like starcheck fids)
    """
    # Only get dwells that have an obsid
    stop = min(DateTime(stop).date, events.obsids.all().reverse()[0].stop)

    logger.info('Getting dwells betweeen {} and {}'.format(start, stop))
    dwells = OrderedDict()
    for dwell in events.dwells.filter(start, stop):
        obsid = dwell.get_obsid()
        if obsid > 40000:
            continue
        if obsid in dwells:
            logger.info('Skipping duplicate obsid {} for dwell {}'.format(obsid, dwell))
            continue
        dwells[obsid] = dwell

    return dwells


def get_fids_telem(dwells, fids_starcheck):
    """
    Get telemetry yag, zag values for each fid in ``fids_starcheck`` at the start
    of the corresponding ``dwells``.

    :returns: dict of tables keyed by obsid (like starcheck fids)
    """

    fids_telem = {}

    for obsid, dwell in dwells.items():
        logger.debug('Get_fids_telem for obsid {}'.format(obsid))
        tstart = dwell.tstart + 300
        if obsid in fids_starcheck:
            fids = fids_starcheck[obsid]
        else:
            logger.info('Skipping obsid {} not in fids_starcheck'.format(obsid))
            continue
        rows = []
        for fid in fids:
            slot = fid['slot']
            try:
                dat = fetch.Msid('aoacyan{}'.format(slot), tstart, tstart + 60)
                yag = np.median(dat.vals)
                dat = fetch.Msid('aoaczan{}'.format(slot), tstart, tstart + 60)
                zag = np.median(dat.vals)
            except Exception as err:
                logger.info('ERROR: {}'.format(err))
            else:
                rows.append((tstart, obsid, slot, yag, zag))

        if len(rows) > 0:
            fids_telem[obsid] = Table(zip(*rows),
                                      names=['tstart', 'obsid', 'slot', 'aoacyan', 'aoaczan'])
        else:
            logger.info('No fids found for obsid {}'.format(obsid))

    return fids_telem


def join_starcheck_telem(fids_starcheck, fids_telem):
    """
    Remake dict of tables into a single table for each structure
    """
    # Stack the dict of tables into a single table
    t_fids_starcheck = table.vstack([fids_starcheck[obsid] for obsid in sorted(fids_starcheck)])
    t_fids_telem = table.vstack([fids_telem[obsid] for obsid in sorted(fids_telem)])

    # Join on obsid and slot columns into a single table
    starcheck_telem = table.join(t_fids_starcheck, t_fids_telem, keys=['obsid', 'slot'])

    # Reject unacquired fids
    ok = starcheck_telem['aoacyan'] > -3276
    return starcheck_telem[ok]


def plot_starcheck_telem(starcheck_telem, savefig=None):
    plt.close(1)
    plt.figure(1, figsize=(6, 4))
    tstart = starcheck_telem['tstart']
    dyag = starcheck_telem['aoacyan'] - starcheck_telem['yang']
    dzag = starcheck_telem['aoaczan'] - starcheck_telem['zang']
    plot_cxctime(tstart, dyag, '.r', label='Yag')
    plot_cxctime(tstart, dzag, '.b', label='Zag')
    plt.ylim(-40, 40)
    plt.grid()
    plt.legend(fontsize='small', numpoints=1)
    plt.ylabel('Offset (arcsec)')
    plt.title('Fid light commanded vs. observed angles')
    x0, x1 = plt.xlim()
    dx = (x1 - x0) * 0.05
    x0, x1 = x0 - dx, x1 + dx
    plt.xlim(x0, x1)
    plt.hlines([-15, 15], x0, x1, colors='g', linestyles='--')
    plt.hlines([-25, 25], x0, x1, colors='r', linestyles='--')
    plt.show()
    if savefig is not None:
        logger.info('Writing {}'.format(savefig))
        plt.savefig(savefig)


def main(args=None):
    opt = get_opt(args)

    start = DateTime(opt.start) - 90 if opt.start is None else DateTime(opt.start)
    stop = DateTime(opt.stop)

    dwells = get_dwells_with_fids(start.date, stop.date)
    fids_starcheck = get_fids_starcheck(dwells)
    fids_telem = get_fids_telem(dwells, fids_starcheck)
    starcheck_telem = join_starcheck_telem(fids_starcheck, fids_telem)
    plot_starcheck_telem(starcheck_telem, savefig=opt.out)
    logger.info('\n'.join(starcheck_telem.pformat(max_lines=-1)))
    return starcheck_telem, fids_starcheck, fids_telem, dwells


if __name__ == '__main__':
    main()
