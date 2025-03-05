#!/usr/bin/env python

import argparse
import json
import os
import re

import matplotlib
import numpy as np
import pyyaks.logger
import tables
from astropy.table import Column, Table
from astropy.time import Time
from kadi import events
from Ska.engarchive import fetch_eng as fetch

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Set up logging
loglevel = pyyaks.logger.INFO
logger = pyyaks.logger.get_logger(
    name="plot_aimpoint", level=loglevel, format="%(asctime)s %(message)s"
)

opt = None  # Global options, set in main

POG = {"ACIS-I": (930.2, 1009.6), "ACIS-S": (200.7, 476.9), "year": 2015.4}

# Aimpoint jumps as seen and defined in fit_aimpoint_drift.ipynb.
AIMPOINT_JUMPS = {
    "2015:006": {"d_dy": 4.4, "d_dz": 1.9},
    "2015:265": {"d_dy": 4.8, "d_dz": 0.5},
    "2016:064": {"d_dy": 2.0, "d_dz": 1.0},
}

acis_arcsec_per_pix = 0.492


def get_opt(args=None):
    parser = argparse.ArgumentParser(
        description="Plot aimpoint drift data from aspect solution files"
    )
    parser.add_argument(
        "--start", default="2000:001", help="Processing start date (default=2000:001)"
    )
    parser.add_argument(
        "--stop", default=Time.now().iso, help="Processing stop date (default=NOW)"
    )
    parser.add_argument(
        "--data-root", default=".", help="Root directory for asol and index files"
    )
    parser.add_argument(
        "--box-duration",
        default=3,
        help="Duration before --stop over which to find aimpoint min/max (months)",
    )

    return parser.parse_args(args)


def get_asol(info=None):
    """
    Get aspect solution DY, DZ, DTHETA values sampled at 1 ksec intervals
    during all science observations.

    :param info: dict of processing information and outputs
    :returns: aspect solution data (Table)
    """
    start = Time(opt.start)
    stop = Time(opt.stop)

    h5_file = os.path.join(opt.data_root, "aimpoint_asol_values.h5")
    logger.info("Reading asol file {}".format(h5_file))
    with tables.open_file(h5_file) as h5:
        asol = Table(h5.root.data[:])
    bad = (asol["dy"] == 0.0) & (asol["dz"] == 0.0)
    asol = asol[~bad]

    year = Column(Time(asol["time"], format="cxcsec").decimalyear, name="year")
    asol.add_column(year, index=0)

    asol.sort("time")

    # Include only points between --start and --stop
    i0, i1 = np.searchsorted(asol["time"], [start.cxcsec, stop.cxcsec])
    asol = asol[i0:i1]

    # Exclude from 10ksec before to 3 days after any normal sun or safe sun intervals.
    normal_suns = events.normal_suns()
    safe_suns = events.safe_suns()
    normal_suns.interval_pad = 10000, 86400 * 3
    safe_suns.interval_pad = 10000, 86400 * 3
    exclude_intervals = (normal_suns | safe_suns).intervals(
        asol["time"][0], asol["time"][-1]
    )
    ok = np.ones(len(asol), dtype=bool)
    for date0, date1 in exclude_intervals:
        logger.info(
            "Excluding asol values from {} to {} (normal/safe sun)".format(date0, date1)
        )
        i0, i1 = np.searchsorted(
            asol["time"], Time([date0, date1], format="yday").cxcsec
        )
        ok[i0:i1] = False
    asol = asol[ok]

    if info is not None:
        info["last_ctime"] = Time(asol["time"][-1], format="cxcsec").datetime.ctime()
        info["last_obsid"] = asol["obsid"][-1]

    for date, jump in AIMPOINT_JUMPS.items():
        jump_date = Time(date)
        if jump_date < stop:
            # Make the mean of the "before" interval match the mean of the "after" interval.
            i0 = np.searchsorted(asol["time"], jump_date.cxcsec)
            asol["dy"][:i0] -= jump["d_dy"] / 20.0
            asol["dz"][:i0] -= jump["d_dz"] / 20.0
            # Capture info about jump
            info.setdefault("aimpoint_jumps", {})[date] = jump
            logger.info("Applying aimpoint jump of {} at {}".format(jump, date))

    return asol


class AsolBinnedStats(object):
    """
    Collect binned statistics from a table of aspect solution offsets.

    This class makes it easy to bin aspect solution dy, dz offsets in
    intervals and then get the percentile values within those bins.

    :param asol: structured ASOL array from update_aimpoint_data.py
    :param bin_days: bin size in days
    :param n_years: use the most recent n_years of data
    """

    def __init__(self, asol, bin_days, n_years=None):
        if n_years is not None:
            iok = np.searchsorted(asol["year"], asol["year"][-1] - n_years)
            asol = asol[iok:]

        t_end = asol[-1]["time"] + 10  # Make sure final bin is just about full
        ibin = (asol["time"] - t_end) // (86400.0 * bin_days)
        self.asol = asol

        for self.det in ("ACIS-S", "ACIS-I"):
            chipx, chipy = self.get_chipx_chipy()
            asol[self.chipx_col] = chipx
            asol[self.chipy_col] = chipy

        self.grouped = asol.group_by(ibin)

    @property
    def chip_col(self):
        return self.det + "_chip"

    @property
    def chipx_col(self):
        return self.det + "_chipx"

    @property
    def chipy_col(self):
        return self.det + "_chipy"

    @property
    def det_title(self):
        return re.sub("-", "", self.det.lower())

    @property
    def ccd(self):
        return "I3" if self.det.lower().endswith("i") else "S3"

    def __getattr__(self, attr):
        """
        If the requested ``attr`` is of the form "p<percentile>_<colname>" for
        an asol column, then compute and return the corresponding percentile
        value for that colname in the grouped asol data structure.
        """
        m = re.match(r"p(\d+)_(\S+)", attr)
        if m:
            perc, col = m.groups()
            _attr = "_" + attr
            if not hasattr(self, _attr):
                rows = []
                for group in self.grouped.groups:
                    isort = np.argsort(group[col])
                    ii = (int(perc) * (len(group) - 1)) // 100
                    rows.append(group[isort[ii]])
                val = Table(rows=rows, names=self.grouped.colnames)
                setattr(self, _attr, val)
            return getattr(self, _attr)
        else:
            return self.__getattribute__(attr)

    def get_chipx_chipy(self):
        if self.det.lower().endswith("i"):
            cx0, cxy, cxz = 971.91, 0.0, -41.74
            cy0, cyy, cyz = 963.07, +41.74, 0.0
        else:
            cx0, cxy, cxz = 252.25, -41.68, 0.0
            cy0, cyy, cyz = 519.95, 0.0, -41.68

        asol = self.asol
        chipx = cx0 + cxy * asol["dy"] + cxz * asol["dz"]
        chipy = cy0 + cyy * asol["dy"] + cyz * asol["dz"]

        return chipx, chipy

    def plot_intra_obs_dy_dz(self):
        # Make a table with the max delta chipx and chipy during each obsid.
        # First get the mins, maxes, means for each obsid by grouping.
        obsids = self.asol.group_by("obsid")["obsid", "year", "dy", "dz"]
        mins = obsids.groups.aggregate(np.minimum)
        maxes = obsids.groups.aggregate(np.maximum)
        means = obsids.groups.aggregate(np.mean)

        t = Table(
            [
                means["year"],
                (maxes["dy"] - mins["dy"]) * 20,
                (maxes["dz"] - mins["dz"]) * 20,
            ],
            names=("year", "dy", "dz"),
        )
        t["ybin"] = np.trunc((t["year"] - t["year"][-1] - 0.0001) * 4.0)

        # Now group the dx and dy vals in 3-month intervals and find the 50th
        # and 90th percentile within each time bin.
        t_bin = t.group_by("ybin")
        outs = {}
        for col in ("dy", "dz"):
            for perc in (50, 90, -5):
                rows = []
                for group in t_bin.groups:
                    i_sort = np.argsort(group[col])
                    if perc < 0:
                        for row in group[i_sort[perc:]]:
                            rows.append(row)
                    else:
                        ii = (int(perc) * (len(group) - 1)) // 100
                        rows.append(group[i_sort[ii]])
                outs[str(perc) + col] = Table(rows=rows, names=t_bin.colnames)

        plt.close(1)
        fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, num=1)
        for col, ax in zip(("dy", "dz"), axes):
            for perc, color in zip((50, 90, -5), ("b", "m", "r")):
                t = outs[str(perc) + col]
                if perc < 0:
                    ax.plot(
                        t["year"],
                        t[col],
                        ".",
                        color=color,
                        markersize=5,
                        label="Outliers",
                    )
                else:
                    ax.plot(
                        t["year"],
                        t[col],
                        color=color,
                        label="{}th percentile".format(perc),
                    )
                    ax.fill_between(t["year"], t[col], color=color, alpha=0.3)
            ax.grid()
            ax.set_ylabel(col.upper() + " (arcsec)")

        ylims = [axes[i].get_ylim()[1] for i in (0, 1)]
        for ii in (0, 1):
            axes[ii].set_ylim(0, max(ylims))
        axes[0].set_title("Intra-observation aimpoint drift (spacecraft coordinates)")
        axes[0].legend(loc="upper left", fontsize="small", title="")
        axes[1].set_xlabel("Year")

        outroot = os.path.join(opt.data_root, "intra_obs_dy_dz")
        logger.info("Writing plot files {}.png".format(outroot))
        fig.patch.set_visible(False)
        plt.savefig(outroot + ".png", facecolor="none")

    def get_chip_x_y_info(self):
        """
        Get information about CHIPX and CHIPY for this detector.

        :param det: detector
        :returns: dict of info
        """
        # Make the N-month bounding box
        asol = self.asol
        iok = np.searchsorted(
            asol["year"], asol["year"][-1] - float(opt.box_duration) / 12
        )
        asol = self.asol[iok:]
        x0, x1 = np.min(asol[self.chipx_col]), np.max(asol[self.chipx_col])
        y0, y1 = np.min(asol[self.chipy_col]), np.max(asol[self.chipy_col])
        ix0, ix1 = np.argmin(asol[self.chipx_col]), np.argmax(asol[self.chipx_col])
        iy0, iy1 = np.argmin(asol[self.chipy_col]), np.argmax(asol[self.chipy_col])
        pogx = POG[self.det][0]
        pogy = POG[self.det][1]

        # Store some information for the web page
        info_det = {}
        info_det["pogx"] = pogx
        info_det["pogy"] = pogy
        info_det["chipx"] = {}
        info_det["chipx"]["min"] = x0
        info_det["chipx"]["mid"] = xmid = (x0 + x1) / 2
        info_det["chipx"]["max"] = x1
        info_det["chipx"]["date_min"] = Time(asol["time"][ix0], format="cxcsec").yday
        info_det["chipx"]["date_max"] = Time(asol["time"][ix1], format="cxcsec").yday
        info_det["chipy"] = {}
        info_det["chipy"]["min"] = y0
        info_det["chipy"]["mid"] = ymid = (y0 + y1) / 2
        info_det["chipy"]["max"] = y1
        info_det["chipy"]["date_min"] = Time(asol["time"][iy0], format="cxcsec").yday
        info_det["chipy"]["date_max"] = Time(asol["time"][iy1], format="cxcsec").yday
        if self.det_title == "aciss":
            info_det["dDY"] = -(pogx - xmid) * acis_arcsec_per_pix
            info_det["dDZ"] = -(pogy - ymid) * acis_arcsec_per_pix
        else:
            info_det["dDY"] = (pogy - ymid) * acis_arcsec_per_pix
            info_det["dDZ"] = -(pogx - xmid) * acis_arcsec_per_pix
        for dd in ("dDY", "dDZ"):
            info_det[dd + "_arcmin"] = info_det[dd] / 60.0

        return info_det

    def plot_chip_x_y(self, info_det):
        """
        Make 3-panel plot showing CHIPX vs. CHIPY, CHIPX vs. time, CHIPY vs. time.

        :param det: detector
        :param info_det: dict of relevant information for plot

        :returns: None
        """
        # Gather plot data from percentile tables
        years = []
        chipxs = []
        chipys = []
        for perc in (0, 10, 50, 90, 100):
            for xy in ("x", "y"):
                dat = getattr(self, "p{}_{}{}".format(perc, self.chip_col, xy))
                years.append(dat["year"])
                chipxs.append(dat[self.chipx_col])
                chipys.append(dat[self.chipy_col])

        year = np.concatenate(years)
        chipx = np.concatenate(chipxs)
        chipy = np.concatenate(chipys)

        # Open three plot axes for CHIPX vs CHIPY, CHIPX vs time and CHIPY vs Time
        plt.close(1)
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2)
        ax2 = plt.subplot2grid((2, 4), (0, 2), colspan=2)
        ax3 = plt.subplot2grid((2, 4), (1, 2), colspan=2, sharex=ax2)

        cm = matplotlib.cm.get_cmap("YlOrRd")

        # Make the N-month bounding box
        asol = self.asol
        iok = np.searchsorted(
            asol["year"], asol["year"][-1] - float(opt.box_duration) / 12
        )
        asol = self.asol[iok:]
        x0, x1 = info_det["chipx"]["min"], info_det["chipx"]["max"]
        y0, y1 = info_det["chipy"]["min"], info_det["chipy"]["max"]
        dx = x1 - x0
        dy = y1 - y0
        year0 = asol["year"][0]
        dyear = asol["year"][-1] - year0
        ax1.add_patch(
            Rectangle((x0, y0), dx, dy, facecolor="#fff8f8", edgecolor="k", zorder=-100)
        )
        ax2.add_patch(
            Rectangle(
                (year0, x0), dyear, dx, facecolor="#fff8f8", edgecolor="k", zorder=-100
            )
        )
        ax3.add_patch(
            Rectangle(
                (year0, y0), dyear, dy, facecolor="#fff8f8", edgecolor="k", zorder=-100
            )
        )
        pogx = POG[self.det][0]
        pogy = POG[self.det][1]

        ax1.plot([pogx], [pogy], "*r", ms=15, zorder=100)
        ax2.plot([POG["year"]], [pogx], "*r", ms=15, zorder=100)
        ax3.plot([POG["year"]], [pogy], "*r", ms=15, zorder=100)

        plot_opt = {
            "c": year,
            "cmap": cm,
            "alpha": 0.8,
            "s": 6.0,
            "linewidths": 0.5,
            "zorder": 10,
        }
        ax1.scatter(chipx, chipy, **plot_opt)
        ax1.set_xlabel("CHIPX")
        ax1.set_ylabel("CHIPY")
        ax1.set_title("{} aimpoint position (CCD {})".format(self.det, self.ccd))
        ax1.set_aspect("equal", "datalim")
        ax1.grid()

        ax2.scatter(year, chipx, **plot_opt)  # points =
        ax2.set_ylabel("CHIPX")
        ax2.yaxis.tick_right()
        ax2.grid()

        ax3.scatter(year, chipy, **plot_opt)  # points =
        ax3.set_xlabel("Year")
        ax3.set_ylabel("CHIPY")
        ax3.yaxis.tick_right()
        ax3.grid()

        outroot = os.path.join(opt.data_root, "chip_x_y_{}".format(self.det_title))
        logger.info("Writing plot files {}.png".format(outroot))
        fig.patch.set_visible(False)
        plt.savefig(outroot + ".png", facecolor="none")


def plot_housing_temperature():
    dat = fetch.Msid("aach1t", "2000:001", stat="daily")
    plt.close(1)
    fig = plt.figure(figsize=(8, 4))
    year = Time(dat.times, format="cxcsec").decimalyear
    plt.plot(year, dat.vals)
    plt.grid()
    plt.xlabel("Year")
    plt.ylabel("Temperature (degF)")
    plt.title("Aspect Camera housing temperature trend")

    outroot = os.path.join(opt.data_root, "aca_housing_temperature")
    logger.info("Writing plot files {}.png".format(outroot))
    fig.patch.set_visible(False)
    plt.savefig(outroot + ".png", facecolor="none")


def make_pure_python(obj):
    """
    Take dict object which can include either dict or numpy scalars or Python scalars, and
    convert to pure Python.
    """
    if isinstance(obj, dict):
        for key, val in obj.items():
            obj[key] = make_pure_python(val)
        return obj
    elif hasattr(obj, "item"):
        return obj.item()
    else:
        return obj


def main():
    global opt
    opt = get_opt()
    info = {
        "date": opt.stop,
        "start": opt.start,
        "stop": opt.stop,
        "box_duration_months": opt.box_duration,
    }

    asol_aimpoint = get_asol(info)

    asol_monthly = AsolBinnedStats(asol_aimpoint, 365.25 / 12)
    for det in ("ACIS-S", "ACIS-I"):
        asol_monthly.det = det
        det_title = asol_monthly.det_title
        info[det_title] = asol_monthly.get_chip_x_y_info()
        asol_monthly.plot_chip_x_y(info[det_title])

    asol_monthly.plot_intra_obs_dy_dz()

    plot_housing_temperature()

    info_file = os.path.join(opt.data_root, "info.json")
    with open(info_file, "w") as fh:
        logger.info("Writing info file {}".format(info_file))
        json.dump(make_pure_python(info), fh, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
