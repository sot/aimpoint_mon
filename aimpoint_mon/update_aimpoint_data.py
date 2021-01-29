#!/usr/bin/env python

import os
import re
import argparse
import tables
from pathlib import Path
import pickle

import numpy as np

from Chandra.Time import DateTime
from astropy.table import Table, Column, vstack
from astropy.time import Time
from mica.archive import asp_l1
from Ska.DBI import DBI
from mica.common import MICA_ARCHIVE
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


def get_asol(obsid, asol_files, dt):
    logger.info('Reading...\n{}'.format('\n'.join(asol_files)))
    asols = [Table.read(asol_file) for asol_file in asol_files]

    # Check to see if the asol files have raw columns ( >= DS 10.8.3)
    has_raws = ['ady' in asol.colnames for asol in asols]
    if np.any(has_raws) and not np.all(has_raws):
        raise ValueError("Some asol files have raw cols and some do not")

    # Reduce to just the columns needed by the tool
    if np.any(has_raws):
        cols = ('time', 'ady', 'adz', 'adtheta')
    else:
        cols = ('time', 'dy', 'dz', 'dtheta')
    asols = [asol[cols] for asol in asols]
    asol = vstack(asols, metadata_conflicts='silent')

    # And rename any raw columns to use the old names
    if np.any(has_raws):
        asol.rename_column('ady', 'dy')
        asol.rename_column('adz', 'dz')
        asol.rename_column('adtheta', 'dtheta')

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
    with tables.open_file(filename, mode='a',
                          filters=tables.Filters(complevel=5, complib='zlib')) as h5:
        try:
            logger.info('Appending {} records to {}'.format(len(asol), filename))
            h5.root.data.append(asol)
        except tables.NoSuchNodeError:
            logger.info('Creating {}'.format(filename))
            h5.createTable(h5.root, 'data', asol, "Aimpoint drift", expectedrows=1e6)
        h5.root.data.flush()


# Set up logging
loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(name='update_aimpoint_data', level=loglevel,
                                  format="%(asctime)s %(message)s")

# Get options
opt = get_opt()
stop = DateTime(opt.stop)
start = stop - 10 if (opt.start is None) else DateTime(opt.start)
logger.info('Processsing from {} to {}'.format(start.date, stop.date))

# Define file names
h5_file = os.path.join(opt.data_root, 'aimpoint_asol_values.h5')


# Get obsids in date range
mica_obspar_db = os.path.join(MICA_ARCHIVE, 'obspar', 'archfiles.db3')
with DBI(dbi='sqlite', server=mica_obspar_db) as db:
    obs = db.fetchall('select obsid, tstart from archfiles where tstart > {}'
                      ' and tstart < {}'
                      .format(start.secs, stop.secs))

# Get unique obsids and then sort by tstart
idx = np.unique(obs['obsid'], return_index=True)[1]
obs = Table(obs[idx])
obs.sort('tstart')
obs['datestart'] = Time(obs['tstart'], format='cxcsec').yday
obs.pprint(max_lines=-1)

# Dict of obsid => list of ASPSOL files. This keeps track of obsids that have
# been processed.
obsid_file = Path(opt.data_root) / 'aimpoint_obsid_index.pkl'
if obsid_file.exists():
    obsid_index = pickle.load(open(obsid_file, 'rb'))
else:
    obsid_index = {}

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

pickle.dump(obsid_index, open(obsid_file, 'wb'))

logger.info('File {} updated'.format(h5_file))
logger.info('File {} updated'.format(obsid_file))

# Write out to FITS
fits_file = re.sub(r'\.h5$', '.fits', h5_file)
dat = Table.read(h5_file, path='data')
dat.meta.clear()
dat.write(fits_file, overwrite=True)
logger.info('File {} updated'.format(fits_file))
