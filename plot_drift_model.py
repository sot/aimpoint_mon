#!/usr/bin/env python
import os
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import Ska.DBI

from Ska.Matplotlib import plot_cxctime
import Ska.engarchive.fetch as fetch
from Chandra.Time import DateTime

plt.rc("axes", labelsize=10)
plt.rc("xtick", labelsize=10)
plt.rc("ytick", labelsize=10)
plt.rc("font", size=10)
plt.rc("legend", fontsize=10)

FIDSET = {'ACIS': [1, 2, 3, 4, 5, 6],
          'HRC-I': [7, 8, 9, 10],
          'HRC-S': [11, 12, 13, 14],
          }


def get_fid_stats(db, det):
    query = ('select id_num, id_string, tstart, ang_y_med, '
             'ang_z_med, sim_z_offset'
             ' FROM fid_stats'
             ' WHERE proc_status IS NULL AND id_string LIKE "{det}%"'
             .format(det=det))
    vals = db.fetchall(query)
    return vals


def plotfids(detstats, det, tstart):
    fidcolor = ' bgrcmkbgrcbgrc'
    fig = plt.figure(1, figsize=(6, 4.0))
    plt.clf()
    ax1 = plt.subplot(1, 1, 1)
    plt.title(det + ' Fid Drift')
    plt.grid()
    plt.ylabel('Y offset (arcsec)')
    plt.grid()
    for fid in FIDSET[det]:
        fidstats = detstats[detstats['id_num'] == fid]
        year = fidstats['tstart'] / 86400. / 365.25 + 1998.0
        normmask = np.logical_and(year > 2002.0, year < 2003.)
        fidstatsnorm = fidstats[normmask]
        if len(fidstatsnorm) < 3:
            print 'Fid %s %d: Not enough readouts for the median' % (det, fid)
            continue
        y0 = np.median(fidstatsnorm['ang_y_med'])
        ok = fidstats['tstart'] > tstart
        if sum(ok) > 0:
            plot_cxctime(fidstats['tstart'][ok],
                         fidstats['ang_y_med'][ok] - y0,
                         'o', markerfacecolor=fidcolor[fid],
                         scaley=False, scalex=False)

    return fig, ax1


def parse_args(args):
    parser = argparse.ArgumentParser(description='Plot drift model')
    parser.add_argument('--data-dir', type=str,
                        default='.',
                        help='Output data directory')
    parser.add_argument('--n-days', type=float,
                        default=100,
                        help='Number of days from present to plot')
    args = parser.parse_args(args)
    return args


def main(args=None):
    args = parse_args(args)
    tstart = DateTime().secs - args.n_days * 86400.0

    db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read')
    detstats = get_fid_stats(db, 'ACIS')
    db.conn.close()
    fig, ax1 = plotfids(detstats, 'ACIS', tstart)

    dat = fetch.MSID('aach1t', tstart, stat='5min')

    yoff = -2.5
    dat0 = 294.6
    adjtemps = -13 - (dat.vals - dat0) * 4 + yoff
    ok = dat.times < DateTime('2011:186').secs
    if sum(ok) > 0:
        plot_cxctime(dat.times[ok], adjtemps[ok], '-r', ax=ax1, fig=fig)
    if sum(~ok) > 0:
        plot_cxctime(dat.times[~ok], adjtemps[~ok], '-r', ax=ax1, fig=fig,
                     alpha=0.3, label='Pre-safemode (2011)')

    adjtemps = -6.5 - (dat.vals - dat0) * 4 + yoff
    ok = dat.times > DateTime('2011:193').secs
    if sum(ok) > 0:
        plot_cxctime(dat.times[ok], adjtemps[ok], '-b', ax=ax1, fig=fig,
                     label='Post-safemode (2011)')
    if sum(~ok) > 0:
        plot_cxctime(dat.times[~ok], adjtemps[~ok], '-b', ax=ax1, fig=fig,
                     alpha=0.3)

    plt.grid()
    plt.legend(loc='upper left')
    plt.tight_layout()

    plt.savefig(os.path.join(args.data_dir, 'drift_model_y.png'))

if __name__ == '__main__':
    main()
