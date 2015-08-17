#!/usr/bin/env python

import re
import os
import argparse
from itertools import izip

import numpy as np
from astropy.table import Table, Column
from astropy.time import Time
import tables
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from bokeh.models import ColumnDataSource
import bokeh.plotting as bp
import mpld3


data_root = '.'

h5_file = os.path.join(data_root, 'aimpoint_asol_values.h5')
POG = {'ACIS-I': (930.2, 1009.6),
       'ACIS-S': (200.7, 476.9),
       'year': 2015.0}


def get_opt(args=None):
    parser = argparse.ArgumentParser(description='Plot aimpoint drift data '
                                     'from aspect solution files')
    parser.add_argument("--data-root",
                        default=".",
                        help="Root directory for asol and index files")
    return parser.parse_args(args)


def get_asol(filename):
    """
    Get aspect solution DY, DZ, DTHETA values sampled at 1 ksec intervals
    during all science observations.
    """
    h5 = tables.openFile(filename)
    asol = Table(h5.root.data[:])
    h5.close()

    bad = (asol['dy'] == 0.0) & (asol['dz'] == 0.0)
    asol = asol[~bad]

    year = Column(Time(asol['time'], format='cxcsec').decimalyear, name='year')
    asol.add_column(year, index=0)

    asol.sort('time')
    return asol


class AsolBinnedStats(object):
    def __init__(self, asol, bin_days):
        t_end = asol[-1]['time'] + 10  # Make sure final bin is just about full
        ibin = (asol['time'] - t_end) // (86400. * bin_days)
        self.asol = asol

        for self.det in ('ACIS-S', 'ACIS-I'):
            chipx, chipy = self.get_chipx_chipy()
            asol[self.chipx_col] = chipx
            asol[self.chipy_col] = chipy

        self.grouped = asol.group_by(ibin)

    @property
    def chip_col(self):
        return self.det + '_chip'

    @property
    def chipx_col(self):
        return self.det + '_chipx'

    @property
    def chipy_col(self):
        return self.det + '_chipy'

    @property
    def det_title(self):
        return re.sub('-', '', self.det.lower())

    @property
    def ccd(self):
        return 'I3' if self.det.lower().endswith('i') else 'S3'

    @property
    def argsort(self):
        if not hasattr(self, '_argsort'):
            self._argsort = self.grouped.groups.aggregate(np.argsort)
        return self._argsort

    def __getattr__(self, attr):
        m = re.match(r'p(\d+)_(\S+)', attr)
        if m:
            perc, col = m.groups()
            _attr = '_' + attr
            if not hasattr(self, _attr):
                rows = []
                for group, isort in izip(self.grouped.groups, self.argsort[col]):
                    ii = (int(perc) * (len(group) - 1)) // 100
                    rows.append(group[isort[ii]])
                val = Table(rows=rows, names=self.grouped.colnames)
                setattr(self, _attr, val)
            return getattr(self, _attr)
        else:
            return self.__getattribute__(attr)

    def get_chipx_chipy(self):
        if self.det.lower().endswith('i'):
            cx0, cxy, cxz = 971.91, 0.0, -41.74
            cy0, cyy, cyz = 963.07, +41.74, 0.0
        else:
            cx0, cxy, cxz = 252.25, -41.68, 0.0
            cy0, cyy, cyz = 963.07, 0.0, -41.68

        asol = self.asol
        chipx = cx0 + cxy * asol['dy'] + cxz * asol['dz']
        chipy = cy0 + cyy * asol['dy'] + cyz * asol['dz']

        return chipx, chipy

    def plot_chip_year(self):
        import bokeh.plotting as bp

        det = self.det

        bp.output_file('chip_year_{}.html'.format(self.det_title),
                       title='{} CHIP vs year'.format(det))
        TOOLS = 'pan,wheel_zoom,box_zoom,crosshair,reset'

        ax1 = bp.figure(tools=TOOLS, toolbar_location='left',
                        plot_width=800, height=300,
                        title='{} aimpoint position (CCD {})'.format(det, self.ccd),
                        y_axis_label='CHIPX')
        ax2 = bp.figure(tools=TOOLS, width=800, height=300,
                        x_range=ax1.x_range,
                        x_axis_label='Year',
                        y_axis_label='CHIPY')

        for ax, xy in izip((ax1, ax2), ('x', 'y')):
            ax.quad(left=self.min['year'],
                    right=self.max['year'],
                    bottom=self.min[self.chip_col + xy],
                    top=self.max[self.chip_col + xy],
                    fill_color='red', fill_alpha=0.5,
                    line_color='black')

            ax.quad(left=self.min['year'],
                    right=self.max['year'],
                    bottom=self.p10[self.chip_col + xy],
                    top=self.p90[self.chip_col + xy],
                    fill_color='blue', fill_alpha=0.5,
                    line_color=None)

        # Put the subplots in a gridplot
        p = bp.gridplot([[ax1],
                         [ax2]])
        bp.save(p)

    def plot_chip_x_y_mpl(self):
        det = self.det

        outfile = 'chip_x_y_{}.html'.format(self.det_title)

        def concat(colname):
            out = [getattr(self, prop)[colname]
                   for prop in ('mean', 'min', 'max')]
            return np.concatenate(out)

        years = []
        chipxs = []
        chipys = []
        for perc in (0, 10, 50, 90, 100):
            for xy in ('x', 'y'):
                dat = getattr(self, 'p{}_{}{}'.format(perc, self.chip_col, xy))
                years.append(dat['year'])
                chipxs.append(dat[self.chipx_col])
                chipys.append(dat[self.chipy_col])

        year = np.concatenate(years)
        chipx = np.concatenate(chipxs)
        chipy = np.concatenate(chipys)

        plt.close(1)
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
        ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)
        ax3 = plt.subplot2grid((2, 4), (1, 2), colspan=2)

        cm = matplotlib.cm.get_cmap('YlOrRd')

        # Make the 6-month bounding box
        asol = self.asol
        iok = np.searchsorted(asol['year'], asol['year'][-1] - 0.5)
        asol = self.asol[iok:]
        x0, x1 = np.min(asol[self.chipx_col]), np.max(asol[self.chipx_col])
        y0, y1 = np.min(asol[self.chipy_col]), np.max(asol[self.chipy_col])
        dx = x1 - x0
        dy = y1 - y0
        ax1.add_patch(Rectangle((x0, y0), dx, dy,
                                facecolor='#fff8f8', edgecolor='k',
                                zorder=-100))
        pogx = POG[self.det][0]
        pogy = POG[self.det][1]
        ax1.plot([pogx], [pogy], '*r', ms=15, zorder=100)
        ax2.plot([POG['year']], [pogx], '*r', ms=15, zorder=100)
        ax3.plot([POG['year']], [pogy], '*r', ms=15, zorder=100)

        opt = dict(c=year, cmap=cm, alpha=0.8, s=6.0, linewidths=0.5, zorder=10)
        points = ax1.scatter(chipx, chipy, **opt)
        ax1.set_xlabel('CHIPX')
        ax1.set_ylabel('CHIPY')
        ax1.set_title('{} aimpoint position (CCD {})'.format(det, self.ccd))
        ax1.set_aspect('equal', 'datalim')
        ax1.grid()

        ax2.scatter(year, chipx, **opt)  # points =
        ax2.set_ylabel('CHIPX')
        ax2.yaxis.tick_right()
        ax2.grid()

        ax3.scatter(year, chipy, **opt)  # points =
        ax3.set_xlabel('Year')
        ax3.set_ylabel('CHIPY')
        ax3.yaxis.tick_right()
        ax3.grid()

        plt.show()
        mpld3.plugins.connect(fig, mpld3.plugins.LinkedBrush(points))

        mpld3.save_html(fig, outfile)

    def plot_chip_x_y(self):
        det = self.det

        bp.output_file('chip_x_y_{}.html'.format(self.det_title),
                       title='CHIP vs year for {}'.format(det))
        TOOLS = 'pan,wheel_zoom,box_zoom,box_select,crosshair,reset'

        year = self.mean['year']
        chipx = self.mean[self.chipx_col]
        chipy = self.mean[self.chipy_col]

        source = ColumnDataSource(data=dict(year=year, chipx=chipx, chipy=chipy))

        width = 500
        ratio = 1.8
        size = 4  # pixels

        ax1 = bp.figure(tools=TOOLS, toolbar_location='left',
                        width=width, height=width,
                        title='{} aimpoint position (CCD {})'.format(det, self.ccd),
                        x_axis_label='CHIPX',
                        y_axis_label='CHIPY')
        ax2 = bp.figure(tools=TOOLS, width=width, height=np.int(width // ratio),
                        x_axis_label='Year',
                        y_axis_label='CHIPX')
        ax3 = bp.figure(tools=TOOLS, width=width, height=np.int(width // ratio),
                        x_axis_label='Year',
                        y_axis_label='CHIPY')

        ax1.circle('chipx', 'chipy', source=source,
                   size=size)

        ax2.circle('year', 'chipx', source=source,
                   size=size)

        ax3.circle('year', 'chipy', source=source,
                   size=size)

        # Put the subplots in a gridplot
        pv = bp.gridplot([[ax2],
                          [ax3]])
        p = bp.gridplot([[ax1], pv])
        bp.save(p)


if 'asol' not in globals():
    asol = get_asol(h5_file)

aw = AsolBinnedStats(asol, 7)
am = AsolBinnedStats(asol, 365.25 / 12)
a3m = AsolBinnedStats(asol, 365.25 / 4)

opt = get_opt(['--data-root=..'])
