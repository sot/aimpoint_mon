import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.table import Table
from matplotlib.colors import to_rgba
from matplotlib.container import ErrorbarContainer

from aimpoint_mon import update_observed_aimpoints
from aimpoint_mon.update_observed_aimpoints import (
    plot_pog_observed_aimpoints,
    select_pog_aimpoints,
    select_recent,
)


@pytest.fixture
def obs_aimpoints() -> Table:
    """
    Small observed-aimpoints table: three in-bounds observations plus one of each
    exclusion case. Obsids 30435 and 32344 are real flight cases with dx ~ 513 arcsec.
    """
    return Table(
        rows=[
            (31955, "ACIS-S", "2026:281:04:50:58.816", 0.84, -1.93, 0.0, 0.0),
            (31686, "ACIS-I", "2026:281:14:35:54.816", -0.75, -2.10, 0.0, 0.0),
            (31506, "HRC-I", "2026:279:20:19:22.816", 15.0, -15.0, 0.0, 0.0),
            (30435, "ACIS-I", "2026:135:22:26:14.816", 513.17, 3.49, 37.8, -39.6),
            (32344, "ACIS-I", "2026:136:12:28:06.816", 1.0, -20.0, 37.8, -39.6),
            (33467, "ACIS-S", "2026:279:14:59:34.816", 0.30, -2.17, 150.0, 0.0),
            (33468, "ACIS-I", "2026:280:18:46:54.816", 0.25, 0.69, 0.0, -150.0),
        ],
        names=[
            "obsid",
            "detector",
            "mean_date",
            "dx",
            "dy",
            "target_offset_y",
            "target_offset_z",
        ],
    )


def test_select_pog_aimpoints_drops_out_of_bounds_from_both_panels(
    obs_aimpoints: Table,
) -> None:
    pog_aimpoints = select_pog_aimpoints(obs_aimpoints)

    # Out of bounds in either dx (30435) or dy (32344) drops the whole observation,
    # as does a target offset over 100 arcsec in Y (33467) or Z (33468). Exactly
    # +/-15 arcsec (31506) is in bounds.
    assert sorted(pog_aimpoints["obsid"]) == [31506, 31686, 31955]


def test_select_recent_keeps_lookback_window_from_newest_observation(
    obs_aimpoints: Table,
) -> None:
    recent = select_recent(obs_aimpoints, lookback_days=30)

    assert sorted(recent["obsid"]) == [31506, 31686, 31955, 33467, 33468]


def test_plot_pog_observed_aimpoints_plots_only_pog_points(
    obs_aimpoints: Table, tmp_path: Path
) -> None:
    output_path = tmp_path / "observed_aimpoints_pog.pdf"

    fig = plot_pog_observed_aimpoints(obs_aimpoints, output_path, lookback_days=180)

    assert output_path.read_bytes().startswith(b"%PDF")
    axes_chipx, axes_chipy = fig.axes
    for axes, expected_offsets in (
        (axes_chipx, [-0.75, 0.84, 15.0]),
        (axes_chipy, [-15.0, -2.10, -1.93]),
    ):
        plotted_offsets = np.concatenate(
            [line.get_ydata() for line in axes.get_lines()]
        )
        assert np.allclose(np.sort(plotted_offsets), expected_offsets)
    assert axes_chipx.get_title() == "Observed aimpoint offsets CHIPX"
    assert axes_chipy.get_title() == "Observed aimpoint offsets CHIPY"
    plt.close(fig)


def test_plot_observed_aimpoints_draws_each_limit_arrow_once_in_detector_color(
    obs_aimpoints: Table, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        update_observed_aimpoints,
        "opt",
        argparse.Namespace(data_root=tmp_path, lookback=180, start=None, stop=None),
    )

    update_observed_aimpoints.plot_observed_aimpoints(obs_aimpoints)

    # The only out-of-range point on each axis is ACIS-I: 30435 has dx = 513 (CHIPX,
    # figure 1) and 32344 has dy = -20 (CHIPY, figure 2). Each is clipped to +/-15
    # and should get one limit arrow, drawn in the ACIS-I color.
    for fig_num, clipped_offset in ((1, 15.0), (2, -15.0)):
        axes = plt.figure(fig_num).axes[0]
        errorbars = [
            container
            for container in axes.containers
            if isinstance(container, ErrorbarContainer)
        ]
        assert len(errorbars) == 1
        data_line, _, (bar_lines,) = errorbars[0].lines
        assert np.allclose(data_line.get_ydata(), [clipped_offset])
        assert to_rgba(data_line.get_color()) == to_rgba("red")
        assert np.allclose(bar_lines.get_color(), [to_rgba("red")])
        plt.close(fig_num)
