#!/usr/bin/env python

import os
import argparse
import functools
from itertools import izip

import numpy as np
from astropy.table import Table, Column
from astropy.time import Time
import matplotlib.pyplot as plt
import tables

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from bokeh.models import ColumnDataSource

data_root = '.'

h5_file = os.path.join(data_root, 'aimpoint_asol_values.h5')


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
    funcs = {'min': np.min,
             'max': np.max,
             'mean': np.mean,
             'p10': functools.partial(np.percentile, q=10),
             'p90': functools.partial(np.percentile, q=90)}
             
    def __init__(self, asol, bin_days):
        t_end = asol[-1]['time'] + 10  # Make sure final bin is just about full
        ibin = (asol['time'] - t_end) // (86400. * bin_days)
        self.asol = asol

        for self.det in ('ACIS-S', 'ACIS-I'):
            chipx, chipy = self.get_chipx_chipy()
            asol[self.chip_col + 'x'] = chipx
            asol[self.chip_col + 'y'] = chipy

        self.grouped = asol.group_by(ibin)

    @property
    def chip_col(self):
        return self.det + '_chip'

    def __getattr__(self, attr):
        if attr in self.funcs:
            _attr = '_' + attr
            if not hasattr(self, _attr):
                val = self.grouped.groups.aggregate(self.funcs[attr])
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

    def plot_chip_vs_time_mpl(self):
        """
        Plot Asol binned stats vs. year
        """
        plt.close(1)
        fig, axes = plt.subplots(nrows=2, sharex=True,
                                 num=1, figsize=(10, 6),
                                 gridspec_kw={'hspace': 0})

        for ax, xy in izip(axes, ('x', 'y')):
            rects = []
            for x0, x1, y0, y1 in izip(self.min['year'], self.max['year'],
                                       self.min[self.chip_col + xy],
                                       self.max[self.chip_col + xy]):
                dx = x1 - x0
                dy = y1 - y0
                r = Rectangle((x0, y0), dx, dy)
                rects.append(r)

            p = PatchCollection(rects, alpha=0.4, facecolors='r',
                                edgecolors='k')
            ax.add_collection(p)
            ax.autoscale()
            ax.grid('on')

        axes[0].set_title('Aimpoint drift for {}'.format(self.det))
        axes[0].set_ylabel('CHIP-X')
        axes[1].set_xlabel('Year')
        axes[1].set_ylabel('CHIP-Y')

        # Turn off overlapping labels
        axes[0].get_yticklabels()[0].set_visible(False)
        axes[1].get_yticklabels()[-1].set_visible(False)

        plt.show()

        return fig, axes

    def plot_chip_year_mission(self):
        import bokeh.plotting as bp

        det = self.det
        ccd = 'I3' if det.lower().endswith('i') else 'S3'

        bp.output_file("chip_year_mission.html",
                       title="CHIP vs year for mission")
        TOOLS = "pan,wheel_zoom,box_zoom,crosshair,reset"

        ax1 = bp.figure(tools=TOOLS, toolbar_location='left',
                        plot_width=800, height=300,
                        title='{} aimpoint position (CCD {})'.format(det, ccd),
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


if 'asol' not in globals():
    asol = get_asol(h5_file)

am = AsolBinnedStats(asol, 365.25 / 12)
a3m = AsolBinnedStats(asol, 365.25 / 4)

opt = get_opt(['--data-root=..'])
