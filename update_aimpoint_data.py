#!/usr/bin/env python

import re
import shelve
import argparse
import os

import tables
import numpy as np
from Chandra.Time import DateTime
from astropy.table import Table, Column, vstack
from astropy.time import Time
from mica.archive import asp_l1
from Ska.DBI import DBI
import pyyaks.logger


def get_opt():
    parser = argparse.ArgumentParser(description='Get aimpoint drift data '
                                     'from aspect solution files')
    parser.add_argument("--data-root",
                        default=".",
                        help="Root directory for asol and index files")
    parser.add_argument("--start",
                        help="Start time for processing (default=stop - 30 days)")
    parser.add_argument("--stop",
                        help="Processing stop date (default=NOW)")
    parser.add_argument("--dt",
                        default=1.0,
                        help="Sample delta time (ksec, default=1.0)")
    return parser.parse_args()


def get_asol(obsid, asol_files, dt, cols=('time', 'dy', 'dz', 'dtheta')):
    logger.info('Reading...\n{}'.format('\n'.join(asol_files)))
    asols = [Table.read(asol_file)[cols] for asol_file in asol_files]
    asol = vstack(asols)

    t0, t1 = asol['time'][[10, -10]]
    n_times = 2 + int((t1 - t0) // (dt * 1000))
    times = np.linspace(t0, t1, n_times)

    idx = np.searchsorted(asol['time'], times)
    asol = asol[idx]
    asol = Table([col.astype(col.dtype.str[1:]) for col in asol.columns.values()])

    asol.add_column(Column([obsid] * len(asol), name='obsid'), index=0)

    return asol


def add_asol_to_h5(filename, asol):
    asol = asol.as_array()
    h5 = tables.openFile(filename, mode='a',
                         filters=tables.Filters(complevel=5, complib='zlib'))
    try:
        logger.info('Appending {} records to {}'.format(len(asol), filename))
        h5.root.data.append(asol)
    except tables.NoSuchNodeError:
        logger.info('Creating {}'.format(filename))
        h5.createTable(h5.root, 'data', asol, "Aimpoint drift", expectedrows=1e6)
    h5.root.data.flush()
    h5.close()


# Set up logging
loglevel = pyyaks.logger.VERBOSE
logger = pyyaks.logger.get_logger(name='get_asol_aimpoint', level=loglevel,
                                  format="%(asctime)s %(message)s")

# Get options
opt = get_opt()
stop = DateTime(opt.stop)
start = stop - 10 if (opt.start is None) else DateTime(opt.start)
logger.info('Processsing from {} to {}'.format(start.date, stop.date))

# Define file names
h5_file = os.path.join(opt.data_root, 'aimpoint_asol_values.h5')
obsid_file = os.path.join(opt.data_root, 'aimpoint_obsid_index.shelve')

# Get obsids in date range
db = DBI(dbi='sqlite', server='/data/aca/archive/obspar/archfiles.db3')
obs = db.fetchall('select obsid, tstart from archfiles where tstart > {}'
                  ' and tstart < {}'
                  .format(start.secs, stop.secs))
db.conn.close()

# Get unique obsids and then sort by tstart
idx = np.unique(obs['obsid'], return_index=True)[1]
obs = Table(obs[idx])
obs.sort('tstart')
obs['datestart'] = Time(obs['tstart'], format='cxcsec').yday
obs.pprint(max_lines=-1)

obsid_index = shelve.open(obsid_file)

# Go through obsids and either process or skip
for obsid in obs['obsid']:
    if str(obsid) in obsid_index:
        logger.info('Skipping obsid {} - already in archive'.format(obsid))
        continue

    logger.info('Processing obsid {}'.format(obsid))
    asol_files = sorted(asp_l1.get_files(obsid=obsid, content='ASPSOL'))
    if not asol_files:
        logger.info('Skipping obsid {} - no asol files'.format(obsid))
        continue

    asol = get_asol(obsid, asol_files, opt.dt)
    add_asol_to_h5(h5_file, asol)
    obsid_index[str(obsid)] = asol_files

obsid_index.close()

logger.info('File {} updated'.format(h5_file))
logger.info('File {} updated'.format(obsid_file))

# Write out to FITS
fits_file = re.sub(r'\.h5$', '.fits', h5_file)
dat = Table.read(h5_file, path='data')
dat.meta.clear()
dat.write(fits_file, overwrite=True)
logger.info('File {} updated'.format(fits_file))
