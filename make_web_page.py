#!/usr/bin/env python

import os
import argparse
import json

from jinja2 import Template
import tables
import pyyaks.logger


def get_opt():
    parser = argparse.ArgumentParser(description='Get aimpoint drift data '
                                     'from aspect solution files')
    parser.add_argument("--data-root",
                        default=".",
                        help="Root directory for asol and index files")
    return parser.parse_args()

# Options
opt = get_opt()

# Set up logging
loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(name='make_web_page', level=loglevel,
                                  format="%(asctime)s %(message)s")
# Files
asol_file = os.path.join(opt.data_root, 'aimpoint_asol_values.h5')
index_template_file = os.path.join(opt.data_root, 'index_template.html')
index_files = {False: os.path.join(opt.data_root, 'index.html'),
               True: os.path.join(opt.data_root, 'index_static.html')}
info_file = os.path.join(opt.data_root, 'info.json')

# Jinja template context
logger.info('Loading info file {}'.format(info_file))
context = json.load(open(info_file, 'r'))

# Get the last record of asol aimpoint values
h5 = tables.open_file(asol_file)
last = h5.root.data[-1]
h5.close()

template = Template(open(index_template_file).read())

for static in (True, False):
    context['static'] = static
    html = template.render(**context)
    index_file = index_files[static]
    logger.info('Writing index file {}'.format(index_file))
    with open(index_file, 'w') as fh:
        fh.write(html)
