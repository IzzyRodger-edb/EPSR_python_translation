"""Generated from fortran_src/plot_pars.inc by tools/gen/gen_common.py -- do not edit.

Common variables are module globals created by ``reset()``; access them as
``plot_pars.<name>``.
"""

import numpy as np  # noqa: F401

from ..fortran.arrays import char_array, dble_array, int8_array, int_array, logical_array, real_array  # noqa: F401
from ..fortran.fmath import f32  # noqa: F401

LOCALS = {}


def reset():
    """(Re)create the COMMON storage of plot_pars.inc."""
    global fname1, fname2, fname1save, fname2save, fnameplot, ext1, ext2, style
    global plabel, title, xlabel, ylabel, zlabel, logx, logy, histogram
    global legend, grid, plottype, plottypeold, nplot, ncolumn, nblock, npcolumn
    global icolumn, ilenfname1, ilenfname2, ilenfnameplot, ilentitle, ilenxlabel, ilenylabel, ilenzlabel
    global ilenlogx, ilenlogy, ilen1save, ilen2save, nplotsave, ncolumnsave, npcolumnsave, xcolumn
    global ilenhistogram, ilenstyle, ilenlegend, ilengrid, icountp, icountf, icountb, ilenext1
    global ilenext2, boffset, boffset2, bspace, xmin, xmax, ymin, yup
    global ymax, ydel, yshift, yfact1, yfact2, xfrac, yfrac
    # common /plot_c/
    fname1 = ''
    fname2 = ''
    fname1save = ''
    fname2save = ''
    fnameplot = ''
    ext1 = ''
    ext2 = ''
    style = ''
    plabel = char_array(100)
    title = ''
    xlabel = ''
    ylabel = ''
    zlabel = ''
    logx = ''
    logy = ''
    histogram = ''
    legend = ''
    grid = ''
    # common /plot_i/
    plottype = 0
    plottypeold = f32(0.0)
    nplot = 0
    ncolumn = 0
    nblock = 0
    npcolumn = 0
    icolumn = int_array(100)
    ilenfname1 = 0
    ilenfname2 = 0
    ilenfnameplot = 0
    ilentitle = 0
    ilenxlabel = 0
    ilenylabel = 0
    ilenzlabel = 0
    ilenlogx = 0
    ilenlogy = 0
    ilen1save = 0
    ilen2save = 0
    nplotsave = 0
    ncolumnsave = 0
    npcolumnsave = 0
    xcolumn = 0
    ilenhistogram = 0
    ilenstyle = 0
    ilenlegend = 0
    ilengrid = 0
    icountp = 0
    icountf = 0
    icountb = 0
    ilenext1 = 0
    ilenext2 = 0
    boffset = 0
    boffset2 = 0
    bspace = 0
    # common /plot_r/
    xmin = f32(0.0)
    xmax = f32(0.0)
    ymin = f32(0.0)
    yup = f32(0.0)
    ymax = f32(0.0)
    ydel = f32(0.0)
    yshift = f32(0.0)
    yfact1 = f32(0.0)
    yfact2 = f32(0.0)
    xfrac = f32(0.0)
    yfrac = f32(0.0)


reset()
