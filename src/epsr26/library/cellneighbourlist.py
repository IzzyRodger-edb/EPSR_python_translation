"""Translation of fortran_src/cellneighbourlist.f.

The routines of ``module cellneighbourlist`` (its data is ``common/cellneighbourlist.py``,
imported here as ``cnl``), followed by the free routines of the file: ``cellcoord``,
``atomcellno``, ``get_date_and_time``, ``get_daysthisyearsofar``, ``get_currentmonth``.
"""

from ..common import cellneighbourlist as cnl  # module cellneighbourlist (data)
from ..common import lattice
from ..common import time as time_common  # common /time/
from ..fortran import Ch, f32, findex, fmt, fstr, idiv, ifix, int_array, nint, powf, r4, real_array, substr, units
from ..fortran import runtime
from ..trace import routine

_F100_init = "(/'initialise_cellneighbourlist> Cell numbers along a,b,c and total',3(1x,i6))"
_F101_init = "(/'initialise_cellneighbourlist> Maximum number of molecules per cell: ',i6)"
_F100_make = "(/'makecellneighbourlist> Found ',i6,' cell neighbours out of ',i6,' total')"

# ***********************************************************************************************


@routine
def initialise_cellneighbourlist(nmolcellin, xadd, yadd, zadd, nmol, rmax):
    """cellneighbourlist.f:initialise_cellneighbourlist -- assign molecules to cells.

    Divides the box into cells holding about nmolcellin molecules each, lists the
    molecules in each cell and the neighbouring cells of the cell at the origin.
    Returns nothing.
    """
    # nmolcell is the (approximate) number of molecules in a cell.
    cnl.nmolcell = nmolcellin
    if cnl.nmolcell > nmol:
        cnl.nmolcell = nmol
    # ncelldivin is the number of cells along the shortest axis. We want the cells to be
    # as close to cubic in each direction so use the Cartesian dimensions to determine cell size
    ncelldiv = nint(powf(f32(nmol) / f32(cnl.nmolcell), r4(1.0) / r4(3.0)))
    if ncelldiv < 1:
        ncelldiv = 1
    cellsize = r4(2.0) * min(lattice.alatx, lattice.blaty, lattice.clatz) / f32(ncelldiv)
    cnl.ncelldiva = 2 * ifix(lattice.alatx / cellsize) + 1
    cnl.ncelldivb = 2 * ifix(lattice.blaty / cellsize) + 1
    cnl.ncelldivc = 2 * ifix(lattice.clatz / cellsize) + 1
    cnl.ncell = cnl.ncelldiva * cnl.ncelldivb * cnl.ncelldivc
    cnl.ncelldivsq = cnl.ncelldiva * cnl.ncelldivb
    units.write(6, fmt(_F100_init, cnl.ncelldiva, cnl.ncelldivb, cnl.ncelldivc, cnl.ncell))
    # Assign the array that will store the number of molecules in each cell.
    if cnl.nmolsincell is not None:
        cnl.nmolsincell = None
    if cnl.molcellno is not None:
        cnl.molcellno = None
    if cnl.iacellneighb is not None:
        cnl.iacellneighb = None
    if cnl.ibcellneighb is not None:
        cnl.ibcellneighb = None
    if cnl.iccellneighb is not None:
        cnl.iccellneighb = None
    cnl.nmolsincell = int_array(cnl.ncell)
    cnl.molcellno = int_array(nmol)
    cnl.iacellneighb = int_array(cnl.ncell)
    cnl.ibcellneighb = int_array(cnl.ncell)
    cnl.iccellneighb = int_array(cnl.ncell)
    for ic in range(1, cnl.ncell + 1):
        cnl.nmolsincell[ic] = 0
    # Assign molecules to cells, storing the number of molecules in each cell
    cnl.maxcellno = 0
    for im in range(1, nmol + 1):
        cnl.molcellno[im] = atomcellno(xadd[im], yadd[im], zadd[im])
        cnl.nmolsincell[cnl.molcellno[im]] = cnl.nmolsincell[cnl.molcellno[im]] + 1
        if cnl.nmolsincell[cnl.molcellno[im]] > cnl.maxcellno:
            cnl.maxcellno = int(cnl.nmolsincell[cnl.molcellno[im]])
    units.write(6, fmt(_F101_init, cnl.maxcellno))
    #     */'Cell number and no. of molecules in cell.')
    # Need to make a list of the molecules in each cell. We assume there will never be more than the current maximum
    # in any one cell
    if cnl.molsincell is not None:
        cnl.molsincell = None
    cnl.molsincell = int_array(cnl.maxcellno, cnl.ncell)
    for ic in range(1, cnl.ncell + 1):
        cnl.nmolsincell[ic] = 0
    for im in range(1, nmol + 1):
        cnl.nmolsincell[cnl.molcellno[im]] = cnl.nmolsincell[cnl.molcellno[im]] + 1
        cnl.molsincell[cnl.nmolsincell[cnl.molcellno[im]], cnl.molcellno[im]] = im
    #      do ic=1,ncell
    #         write(6,*) 'Cell ',ic,(molsincell(i,ic),i=1,nmolsincell(ic))
    #      end do

    # Make the list of neighbouring cells to any one cell at the origin
    makecellneighbourlist(rmax)

    return


# ************************************************************************************************


@routine
def makecellneighbourlist(rmax):
    """cellneighbourlist.f:makecellneighbourlist -- cells within rmax of the origin cell.

    Gets the list of neighbouring cells relative to the one at (0,0,0) which satisfy the
    required distance criterion. Cells which are periodic images of cells already in the
    list are ignored.  Returns nothing.
    """
    # real :: rmax,rmaxsq,rdiffsq !Maximum separation between any atoms in one cell with those in a neighbouring cell
    # real :: xcorner(8),ycorner(8),zcorner(8) !Cartesian coordinates of the corners of each cell relative to the centre of the cell
    xcorner = real_array(8)
    ycorner = real_array(8)
    zcorner = real_array(8)
    # real :: cella,cellb,cellc,cella2,cellb2,cellc2 !Half fractions of cell sides

    cella = r4(1.0) / f32(cnl.ncelldiva)
    cellb = r4(1.0) / f32(cnl.ncelldivb)
    cellc = r4(1.0) / f32(cnl.ncelldivc)
    cella2 = r4(0.5) * cella
    cellb2 = r4(0.5) * cellb
    cellc2 = r4(0.5) * cellc
    # Set up coordinates of corner of cell
    ifound = 0
    for ic in range(-1, 2, 2):
        cx = f32(ic) * cellc2 * lattice.clatc[1]
        cy = f32(ic) * cellc2 * lattice.clatc[2]
        cz = f32(ic) * cellc2 * lattice.clatc[3]
        for ib in range(-1, 2, 2):
            bx = f32(ib) * cellb2 * lattice.blatc[1]
            by = f32(ib) * cellb2 * lattice.blatc[2]
            bz = f32(ib) * cellb2 * lattice.blatc[3]
            for ia in range(-1, 2, 2):
                ax = f32(ia) * cella2 * lattice.alatc[1]
                ay = f32(ia) * cella2 * lattice.alatc[2]
                az = f32(ia) * cella2 * lattice.alatc[3]
                ifound = ifound + 1
                xcorner[ifound] = ax + bx + cx
                ycorner[ifound] = ay + by + cy
                zcorner[ifound] = az + bz + cz
    ncorner = ifound
    #      do i=1,ncorner
    #         write(6,*) i,xcorner(i),ycorner(i),zcorner(i)
    #      end do
    # Now get the list of neighbouring cells, the shortest distance between which is less than rmax
    nrangea = (cnl.ncelldiva - 1) // 2
    nrangeb = (cnl.ncelldivb - 1) // 2
    nrangec = (cnl.ncelldivc - 1) // 2
    ifound = 0
    rmaxsq = rmax * rmax
    for ic in range(-nrangec, nrangec + 1):
        cx = f32(ic) * cellc * lattice.clatc[1]
        cy = f32(ic) * cellc * lattice.clatc[2]
        cz = f32(ic) * cellc * lattice.clatc[3]
        for ib in range(-nrangeb, nrangeb + 1):
            bx = f32(ib) * cellb * lattice.blatc[1]
            by = f32(ib) * cellb * lattice.blatc[2]
            bz = f32(ib) * cellb * lattice.blatc[3]
            for ia in range(-nrangea, nrangea + 1):
                ax = f32(ia) * cella * lattice.alatc[1]
                ay = f32(ia) * cella * lattice.alatc[2]
                az = f32(ia) * cella * lattice.alatc[3]
                xdiffc = ax + bx + cx
                ydiffc = ay + by + cy
                zdiffc = az + bz + cz
                # step through all the possible combinations of corners to determine the shortest distance
                found = False
                if1 = 0
                while if1 < ncorner and not found:
                    if1 = if1 + 1
                    if2 = 0
                    while if2 < ncorner and not found:
                        if2 = if2 + 1
                        xdiff = xdiffc + xcorner[if1] - xcorner[if2]
                        ydiff = ydiffc + ycorner[if1] - ycorner[if2]
                        zdiff = zdiffc + zcorner[if1] - zcorner[if2]
                        rdiffsq = xdiff * xdiff + ydiff * ydiff + zdiff * zdiff
                        # Always include the cell at the origin
                        found = bool(rdiffsq < rmaxsq) or (ia == 0 and ib == 0 and ic == 0)
                if found:
                    #                  write(6,*) ia,ib,ic
                    ifound = ifound + 1
                    cnl.iacellneighb[ifound] = ia
                    cnl.ibcellneighb[ifound] = ib
                    cnl.iccellneighb[ifound] = ic
    cnl.ncellneighb = ifound
    units.write(6, fmt(_F100_make, ifound, cnl.ncell))

    return


@routine
def make_neighbour_list_from_cells(ic, nnbref, nbref):
    """cellneighbourlist.f:make_neighbour_list_from_cells -- molecules in cells near cell ic.

    Fills nbref(1:nnbref) with the molecules in all the neighbouring cells of cell ic.
    Returns nnbref.
    """
    icellref = 0  # set by getcellreference before use
    nnbref = 0
    # Get the coordinates of cell ic
    iaorigin, iborigin, icorigin = findiaibic(ic, 0, 0, 0)
    # Step through all neighbouring cells, and generate the corresponding cell number
    for icref in range(1, cnl.ncellneighb + 1):
        ianeighb = iaorigin + int(cnl.iacellneighb[icref])
        ibneighb = iborigin + int(cnl.ibcellneighb[icref])
        icneighb = icorigin + int(cnl.iccellneighb[icref])
        ianeighb, ibneighb, icneighb, icellref = getcellreference(ianeighb, ibneighb, icneighb, icellref)
        # Get all the molecules in this cell.
        for iref in range(1, int(cnl.nmolsincell[icellref]) + 1):
            nnbref = nnbref + 1
            nbref[nnbref] = cnl.molsincell[iref, icellref]
    return nnbref


@routine
def get_nnbrefmax(nnbrefmax):
    """cellneighbourlist.f:get_nnbrefmax -- largest number of neighbouring molecules.

    To find the maximum number of neighbouring molecules for any cell (including itself).
    Returns nnbrefmax.
    """
    icellref = 0  # set by getcellreference before use
    nnbrefmax = 0
    for ic in range(1, cnl.ncell + 1):
        nnbref = 0
        # Get the coordinates of cell ic
        iaorigin, iborigin, icorigin = findiaibic(ic, 0, 0, 0)
        # Step through all neighbouring cells, and generate the corresponding cell number
        for icref in range(1, cnl.ncellneighb + 1):
            ianeighb = iaorigin + int(cnl.iacellneighb[icref])
            ibneighb = iborigin + int(cnl.ibcellneighb[icref])
            icneighb = icorigin + int(cnl.iccellneighb[icref])
            ianeighb, ibneighb, icneighb, icellref = getcellreference(ianeighb, ibneighb, icneighb, icellref)
            # Get all the molecules in this cell.
            nnbref = nnbref + int(cnl.nmolsincell[icellref])
        if nnbref > nnbrefmax:
            nnbrefmax = nnbref
    return nnbrefmax


@routine
def findiaibic(cellno, ia, ib, ic):
    """cellneighbourlist.f:findiaibic -- cell number -> (a,b,c) cell coordinates.

    Convert cellnumber to the corresponding a,b,c coordinates.  Returns (ia, ib, ic).
    """
    cellno1 = cellno - 1
    ic = idiv(cellno1, cnl.ncelldivsq)
    ib = idiv(cellno1 - ic * cnl.ncelldivsq, cnl.ncelldiva)
    ia = cellno1 - cnl.ncelldiva * (ib + cnl.ncelldivb * ic)

    return ia, ib, ic


@routine
def getcellreference(ia, ib, ic, icellref):
    """cellneighbourlist.f:getcellreference -- cell number of cell (ia,ib,ic).

    Generates a cell reference based on ia,ib,ic, after correcting the latter for
    periodic boundary conditions.  Returns (ia, ib, ic, icellref).
    """
    ia = findperiodiccell(ia, cnl.ncelldiva)
    ib = findperiodiccell(ib, cnl.ncelldivb)
    ic = findperiodiccell(ic, cnl.ncelldivc)
    icellref = 1 + ia + cnl.ncelldiva * (ib + cnl.ncelldivb * ic)

    return ia, ib, ic, icellref


@routine
def findperiodiccell(ix, ncelldiv):
    """cellneighbourlist.f:findperiodiccell -- fold a cell coordinate into 0..ncelldiv-1.

    Gets the periodic cell coordinate based number of cell divisions.  Returns ix.
    """
    while ix < 0:
        ix = ix + ncelldiv
    while ix >= ncelldiv:
        ix = ix - ncelldiv

    return ix


# end module cellneighbourlist

# ************************************************************************************************


@routine
def cellcoord(x, ncelldiv):
    """cellneighbourlist.f:cellcoord -- cell coordinate of fractional coordinate x.

    Returns the appropriate cell coordinate, applying periodic boundary conditions
    (x is replaced by its minimum-image value).  Returns (cellcoord, x).
    """
    x1 = x - f32(nint(x))
    ix = ifix((x1 + r4(0.5)) * f32(ncelldiv))
    while ix < 0:
        ix = ix + ncelldiv
    while ix >= ncelldiv:
        ix = ix - ncelldiv
    x = x1
    return ix, x


# ***********************************************************************************************


@routine
def atomcellno(x, y, z):
    """cellneighbourlist.f:atomcellno -- cell number of Cartesian position (x,y,z).

    Returns atomcellno.
    """
    # From input Cartesian coordinates, gets cell number

    # First convert coordinates to fractional coordinates along a,b and c axes
    fc = z / lattice.clatz
    fb = (y - fc * lattice.claty) / lattice.blaty
    fa = (x - fb * lattice.blatx - fc * lattice.clatx) / lattice.alatx
    ia, fa = cellcoord(fa, cnl.ncelldiva)
    ib, fb = cellcoord(fb, cnl.ncelldivb)
    ic, fc = cellcoord(fc, cnl.ncelldivc)

    result = 1 + ia + cnl.ncelldiva * (ib + cnl.ncelldivb * ic)
    #      write(6,*) ia,ib,ic,atomcellno,ncelldiva,ncelldivb,ncelldivc,ncell

    return result


# ************************************************************************************************

_F701 = "(i2,'-',a3,'-',i4)"
_F702 = "(i2.2,':',i2.2,':',i2.2,'.',i3.3)"
_F703 = "(a,2x,a,1x,a)"
_F704 = "(1x,'Elapsed time (secs): ',i10,'.',i3.3)"


@routine
def get_date_and_time(text, start):
    """cellneighbourlist.f:get_date_and_time -- print the date and time, or the elapsed time.

    Writes text(1:index(text,'#')-1), the time and the date.  start=0 initialises the
    start time (common /time/); otherwise the elapsed time since the previous call is
    written as well.  Returns nothing.
    """
    _date, _time, _zone, values = runtime.date_and_time()
    #      write(6,'(a,8(1x,i5))') 'get_date_and_time> ',(values(i),i=1,8)

    currentyear = values[1]
    month = get_currentmonth(values[2])

    date = fstr(fmt(_F701, values[3], Ch(month, 3), values[1]), 11)
    time = fstr(fmt(_F702, values[5], values[6], values[7], values[8]), 12)
    #      write(6,'(a,1x,a,1x,a,1x,a)') 'get_date_and_time> ',text(1:len_trim(text)),date,time
    #      write(6,'(a,1x,a,1x,a,1x,a)') 'get_date_and_time> ',text(1:index(text,'#')-1),date,time
    units.write(6, fmt(_F703, substr(text, 1, findex(text, "#") - 1), Ch(time, 12), Ch(date, 11)))

    # Generate the total number of days this year so far, not including the current day

    ndays, values[2] = get_daysthisyearsofar(currentyear, values[2])
    daysthisyearsofar = ndays + values[3] - 1
    hoursthisyearsofar = daysthisyearsofar * 24 + values[5]
    minsthisyearsofar = hoursthisyearsofar * 60 + values[6]
    secsthisyearsofar = minsthisyearsofar * 60 + values[7]
    musecs = values[8]
    #      write(6,*) get_daysthisyearsofar(year,values(2)),daysthisyearsofar

    # Calculate elapsed time

    if start == 0:
        # Initialise start time

        time_common.startyear = currentyear
        time_common.startabsolutetime = secsthisyearsofar
        time_common.startmusecs = values[8]

    else:
        # Calculate the elapsed time

        # First calculate the effect if the year number has changed between start and finish
        # We count the number of seconds from the beginning of the start year to the end of the
        # year before the current year

        year = time_common.startyear
        sumdays = 0
        while year < currentyear:
            ndays, _latest = get_daysthisyearsofar(year, 13)
            sumdays = sumdays + ndays
            year = year + 1
        # Convert this to seconds
        sumsecs = sumdays * (24 * 3600)
        # Subtract the previous absolute time and add on the number of secs this year
        sumsecs = sumsecs - time_common.startabsolutetime + secsthisyearsofar
        if time_common.startmusecs > musecs:
            sumsecs = sumsecs - 1
            musecs = 1000 + musecs - time_common.startmusecs
        else:
            musecs = musecs - time_common.startmusecs

        units.write(6, fmt(_F704, sumsecs, musecs))

        # Store the new values of startyear and startabsolute time
        time_common.startyear = currentyear
        time_common.startabsolutetime = secsthisyearsofar
        time_common.startmusecs = values[8]

    return


@routine
def get_daysthisyearsofar(year, latestmonth):
    """cellneighbourlist.f:get_daysthisyearsofar -- days in year before month latestmonth.

    Calculates the number of days in the specified year up to but not including
    latestmonth. To get all months set latestmonth > 12.
    Returns (get_daysthisyearsofar, latestmonth).
    """
    if latestmonth > 12:
        latestmonth = 13
    i = 0
    daysthisyearsofar = 0
    while i < latestmonth - 1:
        dayspermonth = 0
        i = i + 1
        if i == 1:
            dayspermonth = 31
        elif i == 2:
            # Check for leap year
            if (year - 4 * idiv(year, 4)) == 0:
                dayspermonth = 29
            else:
                dayspermonth = 28
        elif i == 3:
            dayspermonth = 31
        elif i == 4:
            dayspermonth = 30
        elif i == 5:
            dayspermonth = 31
        elif i == 6:
            dayspermonth = 30
        elif i == 7:
            dayspermonth = 31
        elif i == 8:
            dayspermonth = 31
        elif i == 9:
            dayspermonth = 30
        elif i == 10:
            dayspermonth = 31
        elif i == 11:
            dayspermonth = 30
        elif i == 12:
            dayspermonth = 31
        daysthisyearsofar = daysthisyearsofar + dayspermonth

    return daysthisyearsofar, latestmonth


@routine
def get_currentmonth(values):
    """cellneighbourlist.f:get_currentmonth -- three-letter name of month ``values``.

    Returns get_currentmonth (character*3).
    """
    # NB: month is not initialised in the Fortran when values is outside 1..12; blank here
    month = ""

    # Get the current month

    if values == 1:
        month = "JAN"
    elif values == 2:
        month = "FEB"
    elif values == 3:
        month = "MAR"
    elif values == 4:
        month = "APR"
    elif values == 5:
        month = "MAY"
    elif values == 6:
        month = "JUN"
    elif values == 7:
        month = "JUL"
    elif values == 8:
        month = "AUG"
    elif values == 9:
        month = "SEP"
    elif values == 10:
        month = "OCT"
    elif values == 11:
        month = "NOV"
    elif values == 12:
        month = "DEC"

    return month
