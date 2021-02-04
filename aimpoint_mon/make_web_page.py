#!/usr/bin/env python

import os
import argparse
import json
from pathlib import Path

from jinja2 import Template
import pyyaks.logger


def get_opt():
    parser = argparse.ArgumentParser(description='Make aimpoint monitor web page')
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


def main():
    # Files
    index_template_file = Path(__file__).parent / 'data' / 'index_template.html'
    index_file = os.path.join(opt.data_root, 'index.html')
    info_file = os.path.join(opt.data_root, 'info.json')

    # Jinja template context
    logger.info('Loading info file {}'.format(info_file))
    context = json.load(open(info_file, 'r'))

    template = Template(open(index_template_file).read())

    context['static'] = True
    html = template.render(**context)
    logger.info('Writing index file {}'.format(index_file))
    with open(index_file, 'w') as fh:
        fh.write(html)


if __name__ == '__main__':
    main()
