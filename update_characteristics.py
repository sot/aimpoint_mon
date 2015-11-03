#!/usr/bin/env python

"""
Update characteristics file.
"""
from __future__ import print_function, division

import re
import os
from datetime import datetime
import argparse
import json
import copy
from difflib import HtmlDiff

import numpy as np
from astropy.time import Time
from astropy.table import Table, Column
from kadi import occweb
from kadi.events import scrape
from bs4 import BeautifulSoup as parse_html
import pyyaks.logger

from calc_si_align import calc_si_align


RE_ODB_SI_ALIGN = re.compile(r'\s* ODB_SI_ALIGN \s* =', re.VERBOSE)
CHARACTERISTICS_URL = ('http://occweb.cfa.harvard.edu/occweb/FOT/configuration/'
                       'documents/Characteristics_Constraints/')

loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(name=__file__, level=loglevel,
                                  format="%(asctime)s %(message)s")

# Global which are set in main
opt = None
version = None


def get_opt(args=None):
    """
    Get runtime options
    """
    parser = argparse.ArgumentParser(description='Create an updated characterstics, '
                                     'place file on lucky, and send notification emails')
    parser.add_argument("--data-root",
                        default=".",
                        help="Root directory for asol and index files")
    parser.add_argument("--process-time",
                        type=str,
                        help="Specify process time (mostly for testing, default=Now)")
    return parser.parse_args(args)


def get_aimpoint_offsets():
    """
    Get most recent aimpoint offset values

    :returns: tuple of dy_acis_i, dz_acis_i, dy_acis_s, dz_acis_s (arcsec)
    """
    info_file = os.path.join(opt.data_root, 'info.json')
    with open(info_file, 'r') as fh:
        info = json.load(fh)

    process_time = Time(opt.process_time) if opt.process_time else Time.now()
    if (process_time - Time(info['date'])).jd > 14:
        logger.info('WARNING: offsets are more than 2 weeks old, last updated {}'
                    .format(info['date']))

    offsets = (info['acisi']['dDY'], info['acisi']['dDZ'],
               info['aciss']['dDY'], info['aciss']['dDZ'])
    logger.info('Read {} updated {} and found offsets {:.2f}, {:.2f}, {:.2f}, {:.2f}'
                .format(info_file, info['date'], *offsets))

    return offsets


def make_updated_characteristics(baseline_file, offsets=(-5, 10, 15, -20)):
    """
    Make an updated characteristics file, starting from path ``baseline_file``
    and updating ODB_SI_ALIGN based on the supplied DY and DZ offsets for
    ACIS-I and ACIS-S.

    The updated file is named using CHARACTERIS_DDMMMYY based on the current date,
    and stored in the characteristics/ directory of opt.data_root.

    :param baseline_file: baseline file path
    :param offsets: tuple of dy_acis_i, dz_acis_i, dy_acis_s, dz_acis_s (arcsec)

    :returns: info_table
    """
    process_time = Time(opt.process_time) if opt.process_time else Time.now()

    # Calculate new SI_ALIGN matrices based on supplied aimpoint offsets
    dy_acis_i, dz_acis_i, dy_acis_s, dz_acis_s = offsets
    si_align_acis_i = calc_si_align(dy_acis_i / 3600, dz_acis_i / 3600.)
    si_align_acis_s = calc_si_align(dy_acis_s / 3600., dz_acis_s / 3600.)

    logger.info('Reading baseline characteristics {}'.format(baseline_file))
    with open(baseline_file, 'r') as fh:
        lines = fh.readlines()

    # Get the starting line number of the ODB_SI_ALIGN definition
    matches = [i for i, line in enumerate(lines) if RE_ODB_SI_ALIGN.match(line)]
    if len(matches) != 1:
        raise ValueError('parsing error, matched {} instances of ODB_SI_ALIGN instead of one'
                         .format(len(matches)))
    start = matches[0]

    comments = ['!Updated via aimpoint_mon/update_characteristics version {}'.format(version),
                '!Run at {}Z'.format(process_time.iso[:16]),
                '!Started with baseline file {}'.format(os.path.basename(baseline_file)),
                '!ACIS-I offsets DY={:.2f} arcsec, DZ={:.2f} arcsec'.format(dy_acis_i, dz_acis_i),
                '!ACIS-S offsets DY={:.2f} arcsec, DZ={:.2f} arcsec'.format(dy_acis_s, dz_acis_s)]

    fmt_first = '       ODB_SI_ALIGN     = {:.6e}, {:.6e}, {:.6e},'
    fmt_other = '                          {:.6e}, {:.6e}, {:.6e},'

    # Make updated lines in characteristics.  Note the use of the transpose (.T)
    # below so that the output is in Row-major (Fortran) order as required by
    # the CHARACTERISTICS.
    for i, vals in enumerate(np.concatenate([si_align_acis_i.T, si_align_acis_s.T])):
        fmt = fmt_first if i == 0 else fmt_other
        out = fmt.format(*vals).upper()
        if i < len(comments):
            out = out + ' ' * max(72 - len(out), 1) + comments[i]

        lines[i + start] = out + os.linesep

    # Write the new file into the characteristics/ directory
    filename = 'CHARACTERIS_{}'.format(process_time.datetime.strftime('%d%b%g').upper())
    updated_file = os.path.join(opt.data_root, 'characteristics', filename)
    logger.info('Writing updated characteristics {}'.format(updated_file))
    with open(updated_file, 'w') as fh:
        fh.writelines(lines)

    # Info fields for summary index table
    info_table = {'date': process_time.yday[:17],
                  'dy_acis_i': dy_acis_i,
                  'dz_acis_i': dz_acis_i,
                  'dy_acis_s': dy_acis_s,
                  'dz_acis_s': dz_acis_s,
                  'baseline_file': os.path.basename(baseline_file),
                  'updated_file': os.path.basename(updated_file),
                  }

    # More detailed info data for JSON file
    info_json = copy.copy(info_table)
    info_json.update({'date_iso': process_time.isot[:19],
                      'dy_dz_units': 'arcsec',
                      'si_align_acis_i': si_align_acis_i.tolist(),
                      'si_align_acis_s': si_align_acis_s.tolist(),
                      'version': version,
                      })

    # Make a new JSON file with udpate details
    with open(updated_file + '.json', 'w') as fh:
        json.dump(info_json, fh, indent=4, sort_keys=True, separators=(',', ': '))

    return info_table


def write_diff_file(updated_file, baseline_file):
    diff = HtmlDiff()

    updated_file = os.path.join(opt.data_root, 'characteristics', updated_file)
    baseline_file = os.path.join(opt.data_root, 'characteristics', baseline_file)

    with open(updated_file, 'r') as fh:
        updated_lines = fh.readlines()
    with open(baseline_file, 'r') as fh:
        baseline_lines = fh.readlines()

    diff_str = diff.make_file(baseline_lines, updated_lines,
                              fromdesc=os.path.basename(baseline_file),
                              todesc=os.path.basename(updated_file),
                              context=True,
                              numlines=5)

    diff_file = updated_file + '_diff.html'
    logger.info('Writing diff to {}'.format(diff_file))
    with open(diff_file, 'w') as fh:
        fh.write(diff_str)


def write_index_file(info_table):
    """
    Write a summary index file in tabular format from ``info_table`` dict, e.g.

        updated_file       baseline_file           date       dy_acis_i dz_acis_i dy_acis_s dz_acis_s
    ------------------- ------------------- ----------------- --------- --------- --------- ---------
    CHARACTERIS_12OCT15 CHARACTERIS_12MAR15 2015:285:01:21:25     -5.00     10.00     15.00    -20.00
    """
    index_file = os.path.join(opt.data_root, 'characteristics', 'index')
    if os.path.exists(index_file):
        index = Table.read(index_file, format='ascii.fixed_width_two_line', guess=False)
        matching = index['updated_file'] == info_table['updated_file']
        index = index[~matching]
        if np.any(matching):
            logger.info('WARNING: replacing existing entry for updated_file={}'
                        .format(info_table['updated_file']))
    else:
        index = Table(names=['updated_file', 'baseline_file', 'date',
                             'dy_acis_i', 'dz_acis_i', 'dy_acis_s', 'dz_acis_s'],
                      dtype=['S19', 'S19', 'S17', 'f', 'f', 'f', 'f'])
    index.add_row(info_table)
    for colname in ('dy_acis_i', 'dz_acis_i', 'dy_acis_s', 'dz_acis_s'):
        index[colname].format = '.2f'

    logger.info('Writing index file {}'.format(index_file))
    with open(index_file, 'w') as fh:
        fh.writelines(line + os.linesep for line in index.pformat())

    # Write an HTML index file as a table
    def self_link(vals):
        """Turn vals into a list of self-referenced links"""
        return ['<a href="{0}">{0}</a>'.format(val) for val in vals]

    # Turn table entries into HTML links to same
    diffs = self_link(x + '_diff.html' for x in index['updated_file'])
    jsons = self_link(x + '.json' for x in index['updated_file'])
    updated_files = self_link(index['updated_file'])
    baseline_files = self_link(index['baseline_file'])
    index.remove_columns(['updated_file', 'baseline_file'])
    index.add_columns([Column(updated_files, name='updated_file'),
                       Column(baseline_files, name='baseline_file'),
                       Column(jsons, name='JSON_info'),
                       Column(diffs, name='diff')],
                      [0, 0, 0, 0])

    index_file = index_file + '.html'
    logger.info('Writing index.html file {}'.format(index_file))
    index.write(index_file, format='ascii.html')

    # Hack: undo the HTML escaping that table write does.
    # TODO: just do the whole thing with jinja template.
    lines = (re.sub(r'&gt;', '>', line) for line in open(index_file, 'r'))
    lines = [re.sub(r'&lt;', '<', line) for line in lines]
    with open(index_file, 'w') as fh:
        fh.writelines(lines)


def get_baseline_characteristics_file():
    """
    Get the most recent (presumed to be baseline) OFLS characteristics file on the OCCweb
    configuration directory.

    :returns: file path (including characteristics/ subdirectory)
    """
    # Get the directory of available files.  This is an HTML doc which consists of single
    # a list of links.
    logger.info('Getting baseline characteristics file')
    occweb.URLS['char_constr'] = '/occweb/FOT/configuration/documents/Characteristics_Constraints/'
    html = occweb.get_url('char_constr')
    html = scrape.cleantext(html)
    page = parse_html(html)

    # Matches CHARACTERIS_DDMMMYY
    RE_CHAR_FILE = re.compile(r'CHARACTERIS_ (?P<date> \d\d [A-Z]{3} \d\d) $', re.VERBOSE)

    # Find every link tag in the document and inspect every one with
    # href=<valid characteristics name>
    dates = []
    filenames = []
    links = page.findAll('a')
    for link in links:
        # If link reference matches regex then include for processing
        filename = link.attrs['href']
        match = RE_CHAR_FILE.match(filename)
        if match:
            date = datetime.strptime(match.group('date'), '%d%b%y')
            dates.append(Time(date).yday)
            filenames.append(filename)
            logger.info('  Valid file: {}'.format(filename))
        else:
            logger.info('  Skipping invalid href: {}'.format(filename))

    # Get the filename for the most recent file
    filename = filenames[np.argmax(dates)]
    logger.info('{} is the most recent characteristics file'.format(filename))

    occweb.URLS[filename] = occweb.URLS['char_constr'] + '/' + filename
    logger.info('Fetching {}'.format(occweb.URLS[filename]))
    html = occweb.get_url(filename)

    outfile = os.path.join(opt.data_root, 'characteristics', filename)
    if not os.path.exists(outfile):
        with open(outfile, 'w') as fs:
            fs.write(html)

    return outfile


def main():
    global opt, version
    opt = get_opt()
    version = open(os.path.join(opt.data_root, 'VERSION'), 'r').read().strip()

    baseline_file = get_baseline_characteristics_file()
    offsets = get_aimpoint_offsets()
    info = make_updated_characteristics(baseline_file, offsets)
    write_index_file(info)
    write_diff_file(info['updated_file'], info['baseline_file'])


if __name__ == '__main__':
    main()
