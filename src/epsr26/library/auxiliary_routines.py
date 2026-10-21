"""Translation of fortran_src/auxiliary_routines.f: small shared helpers.

String helpers (indexlastdot, get_atomicnumber), lattice-angle conversions, EPSR's own
degree trigonometric functions (cosd, sind, acosd, asind, atan2d), 3-vector helpers, the
``ran1`` random number generator (COMMON /ran1com/) with its state I/O, the column data
file reader/writers (r_int, r_int_titles, w_int, w_intr), Euler rotation, smoothing, and
the rotation-matrix / Clebsch-Gordan functions used by the spherical-harmonic tools.
"""

import math

from ..common import chemical_symbols, ran1com, wrcom
from ..fortran import (
    PI4ATAN,
    Ch,
    FortranEOF,
    FortranIOError,
    FortranStop,
    acosf,
    asinf,
    atan2f,
    cosf,
    expf,
    f32,
    findex,
    fmt,
    fstr,
    idiv,
    int_array,
    logf,
    lw,
    r4,
    real_array,
    sinf,
    sqrtf,
    substr,
    units,
    wrap32,
)
from ..trace import routine
from . import parse


@routine
def indexlastdot(text):
    """auxiliary_routines.f:indexlastdot -- position of the last dot in text.

    If there is no dot, the length of text +1 is returned (as if it had a dot)
    [NB: the code actually returns 0 when there is no dot; reproduced].
    Returns indexlastdot.
    """
    ilentext = len(text)  # len_trim(text)
    idot = findex(substr(text, 1, ilentext), ".")
    indexlastdot = idot
    while idot > 0 and indexlastdot < ilentext:
        idot = findex(substr(text, indexlastdot + 1, ilentext), ".")
        indexlastdot = indexlastdot + idot
    return indexlastdot


@routine
def get_atomicnumber(ats):
    """auxiliary_routines.f:get_atomicnumber -- atomic number of chemical symbol ats.

    Returns get_atomicnumber (0 if not found).
    """
    # Search for atomic number of this atom type
    ir = 0
    found = False
    while not found and ir < chemical_symbols.melements:
        ir = ir + 1
        found = ats == chemical_symbols.chsymbol[ir]
    if found:
        get_atomicnumber = ir
    else:
        get_atomicnumber = 0
    return get_atomicnumber


_F601 = "('convertalphbetagamma> Incompatible axes definitions: ',3(1x,f10.5))"


@routine
def convertalphabetagamma(alpha, beta, gam, phib, thetac, phic):
    """auxiliary_routines.f:convertalphabetagamma -- crystallographic alpha, beta, gamma
    to the spherical polar coordinates used in EPSR.

    The a vector is assumed to lie along the x-axis, and the b vector is in the x-y plane.
    phib, thetac, phic are intent(out); on the error returns the values not yet assigned
    keep the caller's values (passed in).
    Returns (phib, thetac, phic).
    """
    pi = PI4ATAN  # 4.0*atan(1.0)
    radconv = pi / r4(180.0)
    degconv = r4(180.0) / pi
    phib = gam
    betar = beta * radconv
    cx = cosf(betar)
    if abs(cx) < r4(1.0e-4):
        thetac = f32(0.0)
        phic = f32(0.0)
    else:
        gamr = gam * radconv
        bx = cosf(gamr)
        by = sinf(gamr)
        alphar = alpha * radconv
        bdotc = cosf(alphar)
        cy = (bdotc - bx * cx) / by
        if abs(cx) > r4(1.0e-4):
            phicr = atan2f(cy, cx)
            phic = phicr * degconv
            cr = sqrtf(cx * cx + cy * cy)
            if abs(cr) <= 1.0:
                thetacr = asinf(cr)
                thetac = degconv * thetacr
            else:
                units.write(6, fmt(_F601, cr, phib, phic))
                return phib, thetac, phic
        else:
            units.write(6, fmt(_F601, alpha, beta, gam))
            return phib, thetac, phic
    #      write(6,*) alpha,beta,gam
    #      write(6,*) phib,thetac,phic
    return phib, thetac, phic


@routine
def convertphibthetacphic(phib, thetac, phic, alpha, beta, gam):
    """auxiliary_routines.f:convertphibthetacphic -- EPSR spherical polar coordinates to
    crystallographic alpha, beta, gamma.

    The a vector is assumed to lie along the x-axis, and the b vector is in the x-y plane.
    Returns (alpha, beta, gam).
    """
    pi = PI4ATAN  # 4.0*atan(1.0)
    radconv = pi / r4(180.0)
    degconv = r4(180.0) / pi
    gam = phib
    phicr = phic * radconv
    thetacr = thetac * radconv
    cxy = sinf(thetacr)
    cx = cxy * cosf(phicr)
    if abs(cx) < r4(1.0e-4):
        alpha = f32(90.0)
        beta = f32(90.0)
    else:
        betar = acosf(cx)
        beta = betar * degconv
        phibr = phib * radconv
        alphar = acosf(cxy * cosf(phibr - phicr))
        alpha = alphar * degconv
    #      write(6,*) phib,thetac,phic
    #      write(6,*) alpha,beta,gam
    return alpha, beta, gam


# parameter (degradconv=3.141592653/180.0), folded in single precision
_DEGRADCONV = r4(3.141592653) / r4(180.0)


@routine
def cosd(deg):
    """auxiliary_routines.f:cosd -- EPSR's own cosine of an angle in degrees.

    Not the gfortran intrinsic (see fmath.cosd_intrinsic).  Returns cosd.
    """
    cosd = cosf(deg * _DEGRADCONV)
    return cosd


@routine
def sind(deg):
    """auxiliary_routines.f:sind -- EPSR's own sine of an angle in degrees.

    Returns sind.
    """
    sind = sinf(deg * _DEGRADCONV)
    return sind


@routine
def acosd(cosphi):
    """auxiliary_routines.f:acosd -- arc cosine in degrees.  Returns acosd."""
    acosd = acosf(cosphi)
    acosd = acosd / _DEGRADCONV
    return acosd


@routine
def asind(sinphi):
    """auxiliary_routines.f:asind -- arc sine in degrees.  Returns asind."""
    asind = asinf(sinphi)
    asind = asind / _DEGRADCONV
    return asind


@routine
def atan2d(sinphi, cosphi):
    """auxiliary_routines.f:atan2d -- atan2 in degrees.  Returns atan2d."""
    atan2d = atan2f(sinphi, cosphi)
    atan2d = atan2d / _DEGRADCONV
    return atan2d


@routine
def perpd(a, b, c):
    """auxiliary_routines.f:perpd -- perpendicular projection of vector c onto the
    direction normal to the plane defined by vectors a and b.

    Returns perpd.
    """
    t = real_array(3)
    #
    # form vector product a X b
    #
    vecpro(a, b, t)
    #
    # normalise it to a unit vector
    #
    amp = sqrtf(scapro(t, t))
    for i in range(1, 4):
        t[i] = t[i] / amp
    #
    # project c onto this normal vector
    #
    perpd = scapro(c, t)
    return perpd


@routine
def vecpro(a, b, c):
    """auxiliary_routines.f:vecpro -- vector product c = a x b (c filled in place)."""
    # performs the vector product of c1 x c2 and puts result in c3
    c[1] = a[2] * b[3] - a[3] * b[2]
    c[2] = a[3] * b[1] - a[1] * b[3]
    c[3] = a[1] * b[2] - a[2] * b[1]


@routine
def scapro(a, b):
    """auxiliary_routines.f:scapro -- dot product a.b.  Returns scapro."""
    # Performs dot product a.b
    scapro = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    return scapro


@routine
def mulvec(b, a):
    """auxiliary_routines.f:mulvec -- multiplies vector b by factor a (in place)."""
    for i in range(1, 4):
        b[i] = b[i] * a


@routine
def divvec(b, a):
    """auxiliary_routines.f:divvec -- divide vector b by a if a .ne. 0.0 (in place)."""
    if a != 0.0:
        for i in range(1, 4):
            b[i] = b[i] / a


@routine
def normvec(a):
    """auxiliary_routines.f:normvec -- creates a unit vector in the direction of a."""
    amp = sqrtf(scapro(a, a))
    if amp > 0.0:
        for i in range(1, 4):
            a[i] = a[i] / amp


@routine
def initran1(nchan):
    """auxiliary_routines.f:initran1 -- read the ran1 state from unit nchan, or else
    initialise it.

    read(nchan,*,iostat=ierr) mcall,idum,iy,(iv(i),i=1,ntab); the items read before an
    error keep their new values, as in Fortran.
    """
    #
    # open file containing random number parameters, or else initialise
    #
    ierr = 0
    try:
        r = units.list_read(nchan)
        ran1com.mcall = r.int(ran1com.mcall)
        ran1com.idum = r.int(ran1com.idum)
        ran1com.iy = r.int(ran1com.iy)
        for i in range(1, ran1com.ntab + 1):
            ran1com.iv[i] = r.int(ran1com.iv[i])
    except FortranEOF:
        ierr = -1
    except FortranIOError as e:
        ierr = e.iostat
    #
    # if the number of calls to ran1 exceeds 1.0E8, the sequence is reset, as this
    # is approximately 1/20th of the period of the ran1 generator
    #
    #           if(mcall.gt.100000000) then
    #                 idum=0
    #                 iy=0
    #                 mcall=0
    #           endif

    # if error or end of file occurs, initialise values

    if ierr != 0:
        ran1com.idum = 0
        ran1com.iy = 0
        ran1com.mcall = 0


@routine
def saveran1(nchan):
    """auxiliary_routines.f:saveran1 -- write the ran1 state to unit nchan.

    Returns nchan (as SIGNATURES.md lists; it is not changed).
    """
    units.write(nchan, lw(ran1com.mcall, ran1com.idum, ran1com.iy,
                          *[ran1com.iv[i] for i in range(1, ran1com.ntab + 1)]))
    return nchan


@routine
def ran1():
    """auxiliary_routines.f:ran1 -- "Minimal" random number generator of Park and Miller
    with Bays-Durham shuffle and added safeguards (Numerical Recipes).

    Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
    values). Call with idum a negative integer to initialize; thereafter, do not alter
    idum between successive deviates in a sequence. RNMX should approximate the largest
    floating value that is less than 1.  State in COMMON /ran1com/.
    Returns ran1 (real*4).
    """
    IA = ran1com.ia
    IM = ran1com.im
    IQ = ran1com.iq
    IR = ran1com.ir
    NTAB = ran1com.ntab
    iv = ran1com.iv
    idum = ran1com.idum
    iy = ran1com.iy
    #
    # Initialize if idum.le.0
    #
    if idum <= 0 or iy == 0:
        #
        # Be sure to prevent idum = 0.  (-idum wraps for idum = -2**31, as in Fortran)
        #
        idum = max(wrap32(-idum), 1)
        #
        # Load the shuffle table (after 8 warm-ups).
        #
        for j in range(NTAB + 8, 0, -1):
            k = idum // IQ  # idum > 0 here
            idum = IA * (idum - k * IQ) - IR * k
            if idum < 0:
                idum = idum + IM
            if j <= NTAB:
                iv[j] = idum
        iy = int(iv[1])
    #
    # Start here when not initializing.
    #
    k = idum // IQ  # idum > 0 here
    #
    # Compute idum=mod(IA*idum,IM) without overflows by Schrage's method.
    #
    idum = IA * (idum - k * IQ) - IR * k
    if idum < 0:
        idum = idum + IM
    #
    # Will be in the range 1:NTAB.
    #
    j = 1 + iy // ran1com.ndiv  # iy >= 0
    iy = int(iv[j])
    #
    # Output previously stored value and refill the shuffle table.
    #
    iv[j] = idum
    #
    # Because users don't expect endpoint values.
    #
    ran1 = min(ran1com.am * f32(iy), ran1com.rnmx)
    ran1com.idum = idum
    ran1com.iy = iy
    ran1com.mcall = wrap32(ran1com.mcall + 1)
    return ran1


_F101_R = "(1x,'Read ',i4,' points from ',a)"


@routine
def r_int(lpt, fname, nfile, xd, yd, ed, nxdim, nfdim, ierr):
    """auxiliary_routines.f:r_int -- read data written by w_int (fname.<ext>01, 02, ...).

    fname (character*256) is changed: it ends with the number of the last file opened.
    Returns (lpt, fname, ierr).
    """
    limit = r4(1.0e-15)
    nt = idiv(nfile + (wrcom.nwrmax - 1), wrcom.nwrmax)
    ilen = indexlastdot(fname) + 1
    icount = 0
    for it in range(1, nt + 1):
        #
        # number of elements to read for each x value
        #
        iwr = nfile - (it - 1) * wrcom.nwrmax
        if iwr > wrcom.nwrmax:
            iwr = wrcom.nwrmax
        wrcom.ext = fstr(fmt("(i2.2)", it), 2)  # write(ext,331) it
        fname = fstr(substr(fname, 1, ilen) + substr(wrcom.ext, 1, 2), 256)
        try:
            units.open(8, fname, status="old")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        if ierr != 0:
            return lpt, fname, ierr
        #
        # check to see if the first line is a comment line
        #
        sometext = fstr(units.read_line(8), 256)  # read(8,100) sometext
        ilenhash = findex(sometext, "#")
        if ilenhash <= 0:
            #
            # hash not found, therefore assume the line is a line of numbers -
            # start the file from the beginning
            #
            units.rewind(8)
        #
        # READ it
        #
        i = 1
        # 200 continue
        while True:
            try:
                r = units.list_read(8)
                wrcom.xo[i] = r.real(wrcom.xo[i])
                for ir in range(1, iwr + 1):
                    wrcom.yo[i, ir] = r.real(wrcom.yo[i, ir])
                    wrcom.eo[i, ir] = r.real(wrcom.eo[i, ir])
            except FortranEOF:
                break  # go to 201
            # Remove absurdly small values
            for ir in range(1, iwr + 1):
                if abs(wrcom.yo[i, ir]) < limit:
                    wrcom.yo[i, ir] = 0.0
                if abs(wrcom.eo[i, ir]) < limit:
                    wrcom.eo[i, ir] = 0.0
            i = i + 1
        # 201
        lpt = i - 1
        units.write(6, fmt(_F101_R, lpt, substr(fname, 1, ilen + 2)))
        #
        # save the data for this file
        #
        for ir in range(1, iwr + 1):
            icount = icount + 1
            for i in range(1, lpt + 1):
                yd[i, icount] = wrcom.yo[i, ir]
                ed[i, icount] = wrcom.eo[i, ir]
        units.close(8)
    #
    # copy the x values
    #
    for i in range(1, lpt + 1):
        xd[i] = wrcom.xo[i]
    return lpt, fname, ierr


@routine
def r_int_titles(fname, nfile, titlex, titled, nfdim, ierr):
    """auxiliary_routines.f:r_int_titles -- read the column titles in the first row of
    the files fname.<ext>01, 02, ... (character*14 titlex, character*29 titled).

    Returns (fname, nfile, titlex, ierr).
    """
    nworddim = 201
    nstsize = 3020
    ncf = int_array(nworddim)
    ncl = int_array(nworddim)
    nwords = 0

    nt = idiv(nfile + (wrcom.nwrmax - 1), wrcom.nwrmax)  # noqa: F841
    ilen = indexlastdot(fname) + 1
    icount = 0
    it = 0
    ierr = 0
    nwr = wrcom.nwrmax
    # Keep reading either until an error occurs (file does not exist),
    # or the number of titles read for a file is less than nwrmax (A maximum of nwrmax
    # titles can be written for each file).
    while ierr == 0 and nwr >= wrcom.nwrmax:
        it = it + 1
        wrcom.ext = fstr(fmt("(i2.2)", it), 2)
        fname = fstr(substr(fname, 1, ilen) + substr(wrcom.ext, 1, 2), 256)
        #         write(6,*) 'r_int_titles> ',fname(1:len_trim(fname))
        try:
            units.open(8, fname, status="old")
            ierr = 0
        except FortranIOError as e:
            ierr = e.iostat
        if ierr == 0:
            #
            # check to see if the first line is a comment line
            #
            sometext = fstr(units.read_line(8), nstsize)
            if findex(sometext, "#") > 0:
                nwords = parse.parse(sometext, nwords, ncf, ncl, 20, nstsize, nworddim)
                if nwords > 1:
                    titlex = fstr(substr(sometext, ncf[2], ncl[2]), 14)
                    iw = 2
                    ir = 0
                    while iw < nwords and ir < wrcom.nwrmax:
                        iw = iw + 1
                        ir = ir + 1
                        wrcom.titleo[ir] = fstr(substr(sometext, ncf[iw], ncl[iw]), 29)
                    nwr = ir
                    ir = 0
                    while ir < nwr and icount < nfdim:
                        ir = ir + 1
                        icount = icount + 1
                        titled[icount] = wrcom.titleo[ir]
                else:
                    ierr = 1
            else:
                ierr = 1
            units.close(8)
    nfile = icount
    units.write(6, lw("r_int_titles> nfile: ", nfile))
    return fname, nfile, titlex, ierr


_F333 = "(a12,i3.3,a9)"
_F332 = "(a1,i3.3,a11)"
_F331 = "(i2.2)"
_F101_W = "(1x,'Writing ',i4,' points to ',a)"


@routine
def w_int(lpt, fname, nfile, xd, yd, ed, titlex, titled, nxdim, nfdim):
    """auxiliary_routines.f:w_int -- write x, (y, error) columns to fname.<ext>01, 02, ...
    (at most nwrmax y columns per file), with a '#' title line.

    fname (character*256) is changed: it ends with the number of the last file written.
    Returns fname.
    """
    nt = idiv(nfile + (wrcom.nwrmax - 1), wrcom.nwrmax)
    ilen = indexlastdot(fname) + 1
    icount = 0
    #
    # copy the x values
    #
    for i in range(1, lpt + 1):
        wrcom.xo[i] = xd[i]
    for it in range(1, nt + 1):
        #
        # number of elements to write for each x value
        #
        iwr = nfile - (it - 1) * wrcom.nwrmax
        if iwr > wrcom.nwrmax:
            iwr = wrcom.nwrmax
        #
        # setup output format specification
        #
        nblock = 2 * iwr + 1
        formattitle = fstr(fmt(_F333, "(1h#,1x,a14,", iwr, "(1x,a29))"), 30)
        formatout = fstr(fmt(_F332, "(", nblock, "(1x,e14.7))"), 30)
        lrec = nblock * 15 + 1  # noqa: F841  (recl= of the open)
        units.write(6, lw(Ch(formatout, 30)))
        wrcom.ext = fstr(fmt(_F331, it), 2)
        fname = fstr(substr(fname, 1, ilen) + substr(wrcom.ext, 1, 2), 256)
        units.write(6, fmt(_F101_W, lpt, substr(fname, 1, ilen + 2)))
        units.open(8, fname, status="unknown", form="formatted")
        #
        # get the data for this file
        #
        for ir in range(1, iwr + 1):
            icount = icount + 1
            wrcom.titleo[ir] = titled[icount]
            for i in range(1, lpt + 1):
                wrcom.yo[i, ir] = yd[i, icount]
                wrcom.eo[i, ir] = ed[i, icount]
        #
        # write it
        #
        units.write(8, fmt(formattitle, Ch(titlex, 14), *[Ch(wrcom.titleo[ir], 29) for ir in range(1, iwr + 1)]))
        xo = wrcom.xo
        yo = wrcom.yo
        eo = wrcom.eo
        for i in range(1, lpt + 1):
            vals = [xo[i]]
            for ir in range(1, iwr + 1):
                vals.append(yo[i, ir])
                vals.append(eo[i, ir])
            units.write(8, fmt(formatout, *vals))
        units.close(8)
    # 100   format(13(1x,e14.7))
    return fname


@routine
def w_intr(lpt1, lpt2, fname, nfile, xd, yd, ed, titlex, titled, nxdim, nfdim):
    """auxiliary_routines.f:w_intr -- as w_int, for the rows lpt1..lpt2 only.

    The title format is built with nblock (not iwr) repeats; only iwr titles are given
    so the output is the same.
    Returns fname.
    """
    nt = idiv(nfile + (wrcom.nwrmax - 1), wrcom.nwrmax)
    ilen = indexlastdot(fname) + 1
    icount = 0
    lpt = lpt2 - lpt1 + 1
    #
    # copy the x values
    #
    for i in range(lpt1, lpt2 + 1):
        wrcom.xo[i] = xd[i]
    for it in range(1, nt + 1):
        #
        # number of elements to write for each x value
        #
        iwr = nfile - (it - 1) * wrcom.nwrmax
        if iwr > wrcom.nwrmax:
            iwr = wrcom.nwrmax
        #
        # setup output format specification
        #
        nblock = 2 * iwr + 1
        formattitle = fstr(fmt(_F333, "(1h#,1x,a14,", nblock, "(1x,a29))"), 30)
        formatout = fstr(fmt(_F332, "(", nblock, "(1x,e14.7))"), 30)
        lrec = nblock * 15 + 1  # noqa: F841  (recl= of the open)
        units.write(6, lw(Ch(formatout, 30)))
        wrcom.ext = fstr(fmt(_F331, it), 2)
        fname = fstr(substr(fname, 1, ilen) + substr(wrcom.ext, 1, 2), 256)
        units.write(6, fmt(_F101_W, lpt, substr(fname, 1, ilen + 2)))
        units.open(8, fname, status="unknown", form="formatted")
        #
        # get the data for this file
        #
        for ir in range(1, iwr + 1):
            icount = icount + 1
            wrcom.titleo[ir] = titled[icount]
            for i in range(lpt1, lpt2 + 1):
                wrcom.yo[i, ir] = yd[i, icount]
                wrcom.eo[i, ir] = ed[i, icount]
        #
        # write it
        #
        units.write(8, fmt(formattitle, Ch(titlex, 14), *[Ch(wrcom.titleo[ir], 29) for ir in range(1, iwr + 1)]))
        xo = wrcom.xo
        yo = wrcom.yo
        eo = wrcom.eo
        for i in range(lpt1, lpt2 + 1):
            vals = [xo[i]]
            for ir in range(1, iwr + 1):
                vals.append(yo[i, ir])
                vals.append(eo[i, ir])
            units.write(8, fmt(formatout, *vals))
        units.close(8)
    # 100   format(13(1x,e14.7))
    return fname


@routine
def rotateeuler(xin, yin, zin, phi, theta, chi, xout, yout, zout):
    """auxiliary_routines.f:rotateeuler -- rotate (xin,yin,zin) by Euler angles
    (phi, theta, chi) in radians.

    Returns (xout, yout, zout).
    """
    #      write(6,*) 'rotateeuler> alpha,beta,gamma: ',phi,theta,chi
    calpha = cosf(phi)
    salpha = sinf(phi)
    cbeta = cosf(theta)
    sbeta = sinf(theta)
    cgamma = cosf(chi)
    sgamma = sinf(chi)
    # Set up the 3x3 rotation matrix
    ax = calpha * cbeta * cgamma - salpha * sgamma
    bx = -(calpha * cbeta * sgamma) - salpha * cgamma
    cx = calpha * sbeta
    ay = salpha * cbeta * cgamma + calpha * sgamma
    by = -(salpha * cbeta * sgamma) + calpha * cgamma
    cy = salpha * sbeta
    az = -(sbeta * cgamma)
    bz = sbeta * sgamma
    cz = cbeta
    xout = ax * xin + bx * yin + cx * zin
    yout = ay * xin + by * yin + cy * zin
    zout = az * xin + bz * yin + cz * zin

    return xout, yout, zout


@routine
def smoof(nqmin, nqmax, nsmoo, diff):
    """auxiliary_routines.f:SMOOf -- 3-point smoothing keeping endpoints fixed
    (diff modified in place)."""
    n1 = nqmin + 1
    n2 = nqmax - 1
    if n2 < n1 or nsmoo == 0:
        return
    for it in range(1, nsmoo + 1):
        sf = diff[nqmin]
        sm = diff[n1]
        sum1 = r4(0.25) * (sf + sm)
        for i in range(n1, n2 + 1):
            j = i + 1
            sn = diff[j]
            sum2 = r4(0.25) * (sm + sn)
            diff[i] = sum1 + sum2
            sm = sn
            sum1 = sum2
            # 10    CONTINUE
        # 20    CONTINUE


@routine
def smood(nqmin, nqmax, nsmoo, diff):
    """auxiliary_routines.f:SMOOd -- 3-point smoothing keeping integrated area constant
    (diff modified in place)."""
    n1 = nqmin
    n2 = nqmax - 1
    if n2 < n1 or nsmoo == 0:
        return
    for it in range(1, nsmoo + 1):
        sm = diff[n1]
        diff[n2 + 1] = diff[n2]
        sum1 = r4(0.5) * sm
        for i in range(n1, n2 + 1):
            j = i + 1
            sn = diff[j]
            sum2 = r4(0.25) * (sm + sn)
            diff[i] = sum1 + sum2
            sm = sn
            sum1 = sum2
            # 10    CONTINUE
        # 20    CONTINUE


@routine
def smoop(nqmin, nqmax, nsmoo, diff):
    """auxiliary_routines.f:SMOOp -- 3-point smoothing holding the points nqmin-1,
    nqmax+1 to zero (diff modified in place)."""
    n1 = nqmin
    n2 = nqmax
    if n2 - 2 < n1 or nsmoo == 0:
        return
    for it in range(1, nsmoo + 1):
        sm = f32(0.0)
        sn = diff[n1]
        sum1 = r4(0.25) * (sm + sn)
        for i in range(n1, n2 + 1):
            sm = sn
            j = i + 1
            if j <= n2:
                sn = diff[j]
            else:
                sn = f32(0.0)
            sum2 = r4(0.25) * (sm + sn)
            diff[i] = sum1 + sum2
            sum1 = sum2


def _dexp(x):
    """dexp of a real*8 (overflow gives +Infinity as in Fortran)."""
    try:
        return math.exp(x)
    except OverflowError:
        return math.inf


def _dlog(x):
    """dlog of a real*8 (log(0) = -Infinity, log(<0) = NaN, as in Fortran)."""
    if x > 0.0:
        return math.log(x)
    if x == 0.0:
        return -math.inf
    return math.nan


@routine
def dlm(l, m, n, costh):
    """auxiliary_routines.f:dlm -- the rotation matrix dlmn(theta) (uses lngamma_mod).

    Returns dlm (real*4).
    """
    from ..common import lngamma_mod  # Fortran module data (package A)

    gamln = lngamma_mod.gamln
    # calculates the rotation matrix dlmn(theta)

    dlm = f32(0.0)
    l2 = 2 * l
    lpm = l + m
    if lpm > l2 or lpm < 0:
        return dlm
    lmm = l - m
    if lmm > l2 or lmm < 0:
        return dlm
    lpn = l + n
    if lpn > l2 or lpn < 0:
        return dlm
    lmn = l - n
    if lmn > l2 or lmn < 0:
        return dlm
    mmn = m - n
    mgamma = lngamma_mod.mgamma
    if lpm + 1 > mgamma or lmm + 1 > mgamma or lpn + 1 > mgamma or lmn + 1 > mgamma:
        units.write(6, lw("dlm> ", l, m, n, lpm + 1, lmm + 1, lpn + 1, lmn + 1))
        raise FortranStop(0)
    faco = 0.5 * (gamln[lpm + 1] + gamln[lmm + 1] + gamln[lpn + 1] + gamln[lmn + 1])
    kmax = min(lpm, lmn) + 1
    kmin = max(0, mmn) + 1
    if costh >= 1.0:
        cth2 = 1.0
        dlcth2 = 0.0
        sth2 = 0.0  # noqa: F841
        dlsth2 = -1.0e30
    elif costh <= -1.0:
        cth2 = 0.0  # noqa: F841
        dlcth2 = -1.0e30
        sth2 = 1.0  # noqa: F841
        dlsth2 = 0.0
    else:
        # (1.+costh)*0.5 and sqrt are single precision; the result is stored in real*8
        cth2 = float(sqrtf((r4(1.0) + costh) * r4(0.5)))
        dlcth2 = _dlog(cth2)
        sth2 = float(sqrtf((r4(1.0) - costh) * r4(0.5)))
        dlsth2 = _dlog(sth2)
    l1 = lpm + lmn
    dldif = dlsth2 - dlcth2  # noqa: F841
    sums = 0.0
    nsign = (-1) ** (kmin - 1)
    for k1 in range(kmin, kmax + 1):
        k = k1 - 1
        k2 = 2 * k
        lpmmk = lpm - k + 1
        lmnmk = lmn - k + 1
        kmmmn = k - mmn + 1
        fac = (faco + dlcth2 * float(l1 - k2) + dlsth2 * float(k2 - mmn) - float(gamln[lpmmk])
               - float(gamln[lmnmk]) - float(gamln[kmmmn]) - float(gamln[k1]))
        sums = sums + float(nsign) * _dexp(fac)
        nsign = nsign * (-1)
    dlm = f32(sums)
    return dlm


@routine
def maxdlm(l, m, n):
    """auxiliary_routines.f:maxdlm -- the maximum factorial required for the rotation
    matrix dlmn.  Returns maxdlm."""
    maxdlm = 0
    l2 = 2 * l
    lpm = l + m
    if lpm > l2 or lpm < 0:
        return maxdlm
    lmm = l - m
    if lmm > l2 or lmm < 0:
        return maxdlm
    lpn = l + n
    if lpn > l2 or lpn < 0:
        return maxdlm
    lmn = l - n
    if lmn > l2 or lmn < 0:
        return maxdlm
    mmn = m - n
    maxdlm = max(lpm + 1, lmm + 1, lpn + 1, lmn + 1)
    kmax = min(lpm, lmn) + 1
    kmin = max(0, mmn) + 1
    l1 = lpm + lmn  # noqa: F841
    for k1 in range(kmin, kmax + 1):
        k = k1 - 1
        k2 = 2 * k  # noqa: F841
        lpmmk = lpm - k + 1
        lmnmk = lmn - k + 1
        kmmmn = k - mmn + 1
        maxdlm = max(maxdlm, lpmmk, lmnmk, kmmmn, k1)
    return maxdlm


@routine
def cgcoef(l1, l2, l, m1, m2):
    """auxiliary_routines.f:cgcoef -- Clebsch-Gordan coefficient for specified indices
    (uses lngamma_mod).  Returns cgcoef (real*4)."""
    from ..common import lngamma_mod  # Fortran module data (package A)

    gamln = lngamma_mod.gamln

    cgcoef = f32(0.0)
    m = m1 + m2
    # check for closure relation on l's
    ia = l1 + l2 - l + 1
    if ia <= 0:
        return cgcoef
    ib = l + l1 - l2 + 1
    if ib <= 0:
        return cgcoef
    ic = l + l2 - l1 + 1
    if ic <= 0:
        return cgcoef
    if abs(m) > l:
        return cgcoef
    aid = f32(2 * l + 1)
    ie = l1 + l2 + l + 2
    # real e1 = real*8 expression
    e1 = f32(float(gamln[ia]) + float(gamln[ib]) + float(gamln[ic]) - float(gamln[ie]))
    ja = l1 + m1 + 1
    jb = l1 - m1 + 1
    jc = l2 + m2 + 1
    jd = l2 - m2 + 1
    je = l + m + 1
    jf = l - m + 1
    e1 = f32(float(e1) + float(gamln[ja]) + float(gamln[jb]) + float(gamln[jc]) + float(gamln[jd])
             + float(gamln[je]) + float(gamln[jf]))
    e1 = r4(0.5) * (e1 + logf(aid))
    kt1 = l - l2 + m1
    kt2 = l - l1 - m2
    kmin = max(0, -kt1, -kt2) + 1
    kt3 = l1 + l2 - l + 1
    kt4 = l1 - m1 + 1
    kt5 = l2 + m2 + 1
    kmax = min(kt3, kt4, kt5)
    phase = (-1) ** (kmin - 1)
    sums = f32(0.0)
    for k1 in range(kmin, kmax + 1):
        k = k1 - 1
        e2 = f32(float(gamln[k1]) + float(gamln[kt1 + k1]) + float(gamln[kt2 + k1]) + float(gamln[kt3 - k])
                 + float(gamln[kt4 - k]) + float(gamln[kt5 - k]))
        sums = sums + f32(phase) * expf(e1 - e2)
        phase = -phase
    cgcoef = sums
    return cgcoef


@routine
def maxcgcoef(l1, l2, l, m1, m2):
    """auxiliary_routines.f:maxcgcoef -- the maximum factorial needed for the specified
    Clebsch-Gordan coefficient.  Returns maxcgcoef."""
    #      use lngamma_mod

    maxcgcoef = 0
    m = m1 + m2
    # check for closure relation on l's
    ia = l1 + l2 - l + 1
    if ia <= 0:
        return maxcgcoef
    ib = l + l1 - l2 + 1
    if ib <= 0:
        return maxcgcoef
    ic = l + l2 - l1 + 1
    if ic <= 0:
        return maxcgcoef
    if abs(m) > l:
        return maxcgcoef
    aid = f32(2 * l + 1)  # noqa: F841
    ie = l1 + l2 + l + 2
    maxcgcoef = max(maxcgcoef, ia, ib, ic, ie)
    #     e1=gamln(ia)+gamln(ib)+gamln(ic)-gamln(ie)
    ja = l1 + m1 + 1  # noqa: F841
    jb = l1 - m1 + 1
    jc = l2 + m2 + 1
    jd = l2 - m2 + 1
    je = l + m + 1
    jf = l - m + 1
    maxcgcoef = max(maxcgcoef, jb, jc, jd, je, jf)
    #      e1=e1+gamln(ja)+gamln(jb)+gamln(jc)+gamln(jd)+gamln(je)+gamln(jf)
    #      e1=0.5*(e1+log(aid))
    kt1 = l - l2 + m1
    kt2 = l - l1 - m2
    kmin = max(0, -kt1, -kt2) + 1
    kt3 = l1 + l2 - l + 1
    kt4 = l1 - m1 + 1
    kt5 = l2 + m2 + 1
    kmax = min(kt3, kt4, kt5)
    for k1 in range(kmin, kmax + 1):
        k = k1 - 1
        maxcgcoef = max(maxcgcoef, k1, (kt1 + k1), (kt2 + k1), (kt3 - k), (kt4 - k), (kt5 - k))
        #            e2=gamln(k1)+gamln(kt1+k1)+gamln(kt2+k1)+gamln(kt3-k)+gamln(kt4-k)+gamln(kt5-k)
    return maxcgcoef
