"""Translation of fortran_src/mapgrentropy.f: ``program standalone`` and
``subroutine run_mapgrentropy`` -- estimate the radial, spatial and total (orientational)
two-body excess entropy from the spherical-harmonic coefficients written by sharm/sdf, by
Monte Carlo sampling of random relative orientations.

``setup_symbols_mapgrentropy`` is in setup_symbols_generic.f
(library/setup_symbols_generic.py, imported inside the routine: it is translated by
another work package).  smoosharm and sumoverm2 are in plotsharm_routines.f.

``external cosd,sind,...`` is declared but no degree function is called here; sumoverm2
(plotsharm_routines.f) uses gfortran's intrinsics.
"""

from ..common import commands, symbols
from ..fortran import PI4ATAN, Ch, FortranStop, acosf, char_array, dble_array, f32, findex, fmt, fstr, idiv
from ..fortran import int_array, isign, len_trim, logf, lw, r4, read_internal, real_array, sqrtf, substr, units
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import lngamma_mod
from ..library import parse
from ..library import plotsharm_routines
from ..library import read_and_save_files
from ..trace import routine

# parameter(mexpan=10000,mval=10,mpixel=80,mpixel2=(2*mpixel+1),mpixel3=mpixel2**3)
mexpan = 10000
mval = 10
mpixel = 80
mpixel2 = 2 * mpixel + 1
mpixel3 = mpixel2 ** 3

_F92 = "(1x,10(1x,i5))"


def main():
    """mapgrentropy.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_mapgrentropy()
    raise FortranStop(0)


def _m1pow(n):
    """(-1)**n for an integer n (integer result; also right for negative n)."""
    return 1 if n % 2 == 0 else -1


def _read_index(temptitle, ifirst, lentemptitle, ch):
    """index(temptitle(ifirst:lentemptitle),ch)"""
    return findex(substr(temptitle, ifirst, lentemptitle), ch)


@routine
def run_mapgrentropy():
    """mapgrentropy.f:run_mapgrentropy -- entropy estimate from g(r,w), written to
    <input>.out, <input>.s01 and negativegr.txt; the ran1 state is kept in ran1.dat."""
    from ..library import setup_symbols_generic

    # use lngamma_mod
    # real datac(500,mexpan),datas(500,mexpan),rd(500),grw(500,5),grwerr(500,5)
    # 1,zcoeffc(mexpan),zcoeffs(mexpan),zcoeffspc(mexpan),zcoeffsps(mexpan)
    # *,sumc,sums,spsumc,spsums
    datac = real_array(500, mexpan)
    datas = real_array(500, mexpan)
    rd = real_array(500)
    grw = real_array(500, 5)
    grwerr = real_array(500, 5)
    zcoeffc = real_array(mexpan)
    zcoeffs = real_array(mexpan)
    zcoeffspc = real_array(mexpan)
    zcoeffsps = real_array(mexpan)
    # dimension lvl(3,mval),nl(3),ml(3)
    # 1,th(50),ph(50),xl(2),yl(2),temprd(500),tempsp(500),tempto(500)
    # 2,out(0:mpixel3),xx(201),yy(201),l1ex(mexpan)
    # 3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
    # *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
    # *,tempc(500),temps(500)
    # (lvl, nl, ml, th, ph, xl, yl, out, xx, yy, rdcoor, thcoor, phcoor, radbal, rgbbal are not used)
    temprd = real_array(500)
    tempsp = real_array(500)
    tempto = real_array(500)
    l1ex = int_array(mexpan)
    l2ex = int_array(mexpan)
    lex = int_array(mexpan)
    n1ex = int_array(mexpan)
    n2ex = int_array(mexpan)
    tempc = real_array(500)
    temps = real_array(500)
    # character*256 infile,fname,fnameout
    # character*1 ext
    # character*14 titler
    # character*29 titledat(mexpan)
    titledat = char_array(mexpan)
    # character(len=14)                      :: titlex
    # character(len=29), dimension(mexpan)   :: titled
    # character(len=29)                      :: temptitle
    # character(len=2)                       :: text
    # logical                                :: forcemirrorsym,test
    titlex = ""
    titled = char_array(mexpan)
    # external cosd,sind,acosd,asind,atan2d
    # real pi,term,temp,tsumrd,tsumsp,tsumto
    # real molndens,molmass,moltemp,mollambda,molrot,molvib
    # real kineticterm,rotterm,vibterm,transterm
    # integer idot,nfreq,ident,indentsign,leninfile,lenfname
    # integer nsmoord,nsmoosp,nsmooor,nrminsmoo

    # real(kind=8) entropyrd,entropysp,entropyto,sumrd,sumsp,sumto,grw8(500,3),grwerr8(500,3)
    grw8 = dble_array(500, 3)
    grwerr8 = dble_array(500, 3)
    # real minrd,minsp,minto,minrdr,minspr,mintor
    # real maxrd,maxsp,maxto,maxrdr,maxspr,maxtor
    # real thetalspmin,philspmin,phimspmin,thetamspmin,chimspmin
    # real thetaltomin,philtomin,phimtomin,thetamtomin,chimtomin
    # real spnegative(500),tonegative(500),wtnegative !Store for negative values of g(r) as a function of r.
    # This will be used to correct values at same distance which are above 1.0 to compensate.
    spnegative = real_array(500)
    tonegative = real_array(500)

    #
    # setup the mapgrentropy symbols
    #
    setup_symbols_generic.setup_symbols_mapgrentropy()
    #
    # now setup the input file
    #
    nin = 11
    nout = 6
    # NB: ierrout is not initialised in the Fortran; setup_input_file sets it
    ierrout = 0
    nout, ierrout = file_and_directory_commands.setup_input_file(nin, nout, ierrout)
    #
    # if an error occurred on setup, exit now
    #
    if ierrout != 0:
        return
    #
    # read the input file
    #
    nout, ierrout = read_and_save_files.read_input_file(nin, nout, ierrout)
    #
    # if an error occurred on read, exit now
    #
    if ierrout != 0:
        return
    pi = PI4ATAN  # 4.0*atan(1.0)
    nxdim = 500
    nfdim = mexpan
    #
    # initialise log(factorial)
    #
    lngamma_mod.fctrln(250)
    units.open(10, "ran1.dat", status="unknown")
    auxiliary_routines.initran1(10)
    units.close(10)
    idot = auxiliary_routines.indexlastdot(symbols.fnameinp)
    if idot > 0:
        fnameout = fstr(substr(symbols.fnameinp, 1, idot) + "out", 256)
    else:
        fnameout = fstr(substr(symbols.fnameinp, 1, len_trim(symbols.fnameinp)) + ".out", 256)
    units.write(6, fmt("(a,1x,a)", substr(symbols.fnameinp, 1, len_trim(symbols.fnameinp)),
                       substr(fnameout, 1, len_trim(fnameout))))
    for iref in range(1, symbols.nsymbol + 1):
        units.write(6, fmt("(a,1x,a)", Ch(symbols.symbol[iref], 11), symbols.symbol_value[iref]))
    iref = 1
    infile = symbols.symbol_value[iref]  # Name of file containing SH coefficients
    # convert to full name
    leninfile = len_trim(infile)
    fname = ""
    lenfname = 0
    fname, lenfname = file_and_directory_commands.get_full_filename(infile, leninfile, fname, lenfname)
    #      write(6,97)
    # 97    format(1x,'Type no. of expansion coefficients')
    iref = iref + 1
    nexpan = 0
    nexpan = read_internal(symbols.symbol_value[iref]).int(nexpan)
    #
    # clean out data arrays
    #
    for k in range(1, nexpan + 1):
        datac[1:501, k] = 0.0
        datas[1:501, k] = 0.0
    #      write(6,96)
    # 96    format(1x,'Are molecules identical (ident=0 if so, else 1)')
    # A value of ident .lt. 0 will be used to signify that we want molecule 2 to be
    # at the origin
    iref = iref + 1
    ident = 0
    ident = read_internal(symbols.symbol_value[iref]).int(ident)
    #      read(10,*) ident
    # save the sign of ident and keep separate from the absolute value
    if abs(ident) > 0:
        identsign = isign(1, ident)
        ident = abs(ident)
    else:
        identsign = 1
        ident = 0
    kk0 = 1  # noqa: F841
    # Get the titles from the input file
    # NB: ierr is not initialised in the Fortran; r_int_titles sets it
    ierr = 0
    fname, nexpan, titlex, ierr = auxiliary_routines.r_int_titles(fname, nexpan, titlex, titled, nfdim, ierr)
    for k in range(1, nexpan + 1):
        temptitle = fstr(titled[k], 29)
        lentemptitle = len_trim(temptitle)
        #         write(6,*) k,temptitle
        ifirst = 1
        ifirst = _read_index(temptitle, ifirst, lentemptitle, "(") + 1
        ilast = ifirst + _read_index(temptitle, ifirst, lentemptitle, ",") - 2
        #         write(6,*) ifirst,ilast
        l1ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(l1ex[k])
        ifirst = ilast + 2
        ilast = ifirst + _read_index(temptitle, ifirst, lentemptitle, ",") - 2
        l2ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(l2ex[k])
        ifirst = ilast + 2
        ilast = ifirst + _read_index(temptitle, ifirst, lentemptitle, ";") - 2
        lex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(lex[k])
        ifirst = ilast + 2
        ilast = ifirst + _read_index(temptitle, ifirst, lentemptitle, ",") - 2
        n1ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(n1ex[k])
        ifirst = ilast + 2
        ilast = ifirst + _read_index(temptitle, ifirst, lentemptitle, ")") - 2
        n2ex[k] = read_internal(substr(temptitle, ifirst, ilast)).int(n2ex[k])
        units.write(6, fmt(_F92, k, l1ex[k], l2ex[k], lex[k], n1ex[k], n2ex[k]))
    #
    # read in data to be fitted
    #
    # NB: nr is not initialised in the Fortran; r_int sets it
    nr = 0
    nr, fname, ierr = auxiliary_routines.r_int(nr, fname, nexpan, rd, datac, datas, nxdim, nfdim, ierr)
    iref = iref + 1
    commands.sometext = fstr(symbols.symbol_value[iref], 3020)
    commands.nwords = parse.parse(commands.sometext, commands.nwords, commands.ncf, commands.ncl, 20,
                                  commands.nstsize, commands.nworddim)
    rminsmoo = f32(0.0)
    nsmoord = 0
    nsmoosp = 0
    nsmooor = 0
    if commands.nwords >= 4:
        r = read_internal(commands.sometext)
        rminsmoo = r.real(rminsmoo)
        nsmoord = r.int(nsmoord)
        nsmoosp = r.int(nsmoosp)
        nsmooor = r.int(nsmooor)
    elif commands.nwords >= 3:
        r = read_internal(commands.sometext)
        nsmoord = r.int(nsmoord)
        nsmoosp = r.int(nsmoosp)
        nsmooor = r.int(nsmooor)
        rminsmoo = f32(0.0)
    else:
        nsmoord = read_internal(commands.sometext).int(nsmoord)
        rminsmoo = f32(0.0)
        nsmoosp = nsmoord
        nsmooor = nsmoord
    # Determine the minimum radius value at which to start smoothing.
    nrminsmoo = 1
    while rd[nrminsmoo] < rminsmoo:
        nrminsmoo = nrminsmoo + 1
    # Perform the smoothing on the coefficients as required
    for k in range(1, nexpan + 1):
        #         write(6,92) l1ex(k),l2ex(k),lex(k),n1ex(k),n2ex(k)
        nsmoo = nsmoord
        if l2ex[k] == 0 and l1ex[k] != 0:
            nsmoo = nsmoosp
        elif l2ex[k] != 0 and l1ex[k] != 0:
            nsmoo = nsmooor
        if nsmoo != 0:
            for i in range(1, nr + 1):
                fact = rd[i] * rd[i]
                tempc[i] = datac[i, k] * fact
                temps[i] = datas[i, k] * fact
            plotsharm_routines.smoosharm(nrminsmoo, nr, nsmoo, tempc)
            plotsharm_routines.smoosharm(nrminsmoo, nr, nsmoo, temps)
            for i in range(1, nr + 1):
                fact = rd[i] * rd[i]
                datac[i, k] = tempc[i] / fact
                datas[i, k] = temps[i] / fact
    #      write(6,'(a)') 'Apply mirror symmetry (y or no)'
    #      read(10,'(a)') text
    iref = iref + 1
    text = fstr(symbols.symbol_value[iref], 2)  # read(symbol_value(iref),'(a)') text
    forcemirrorsym = text == "y"
    #      write(6,105)
    # 105   format(1x,'Specify number of trials to calculate entropy for, and frequency to write')
    iref = iref + 1
    nvary = 0
    nvary = read_internal(symbols.symbol_value[iref]).int(nvary)
    iref = iref + 1
    nfreq = 0
    nfreq = read_internal(symbols.symbol_value[iref]).int(nfreq)
    iref = iref + 1
    molndens = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    molmass = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    moltemp = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    molrot = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    iref = iref + 1
    molvib = read_internal(symbols.symbol_value[iref]).real(f32(0.0))
    gasconstant = r4(8.314472)
    # De Broglie wavelength [A]
    mollambda = r4(66.2606876) / sqrtf(f32(2) * pi * molmass * r4(1.660539) * r4(1.3806503) * moltemp)
    transterm = gasconstant * (r4(2.5) - logf(molndens * (mollambda * (mollambda * mollambda))))
    rotterm = gasconstant * (r4(1.5) + logf(sqrtf(pi * (moltemp * moltemp * moltemp) / molrot)))
    vibterm = molvib
    kineticterm = transterm + rotterm + vibterm
    #      read(10,*) nvary,nfreq
    # 121   format(a50)
    # calculate g(r,w) at specified angular values
    entropyrd = float(kineticterm)  # Radial entropy
    entropysp = float(kineticterm)  # Spatial entropy
    entropyto = float(kineticterm)  # total entropy
    for it in range(1, 4):
        for i in range(1, nr + 1):
            grw8[i, it] = float(kineticterm)
            grwerr8[i, it] = float(kineticterm)
    # zero negative value accumulators
    for i in range(1, nr + 1):
        spnegative[i] = 0.0
        tonegative[i] = 0.0
    fact = r4(2.0) * pi * gasconstant * molndens * (rd[3] - rd[2])
    units.open(10, fnameout, status="unknown")
    units.write(10, fmt('("Temperature [K]: ",e12.5)', moltemp))
    units.write(10, fmt('("Number density [mols/A**3]: ",e12.5)', molndens))
    units.write(10, fmt('("De Broglie wavelength: ",e12.5)', mollambda))
    units.write(10, fmt('("Perfect gas entropy, trans. term: ",e12.5)', transterm))
    units.write(10, fmt('("Perfect gas entropy, rot. term: ",e12.5)', rotterm))
    units.write(10, fmt('("Perfect gas entropy,vib. term: ",e12.5)', vibterm))
    units.write(10, fmt('("Perfect gas entropy, total: ",e12.5)', kineticterm))
    minrd = f32(0.0)
    minsp = f32(0.0)
    minto = f32(0.0)
    maxrd = f32(0.0)
    maxsp = f32(0.0)
    maxto = f32(0.0)
    # NB: the radius and angle values of the minima/maxima are not initialised in the
    # Fortran; they are only set when a minimum/maximum is found.  0 reproduces the
    # reference (Rust mapgrentropy.rs).
    minrdr = minspr = mintor = f32(0.0)
    maxrdr = maxspr = maxtor = f32(0.0)
    thetalspmin = philspmin = phimspmin = thetamspmin = chimspmin = f32(0.0)
    thetaltomin = philtomin = phimtomin = thetamtomin = chimtomin = f32(0.0)
    thetalspmax = philspmax = phimspmax = thetamspmax = chimspmax = f32(0.0)
    thetaltomax = philtomax = phimtomax = thetamtomax = chimtomax = f32(0.0)
    l1only = 1
    l2only = 1
    lonly = 1
    n1only = 1
    n2only = 1
    m2fix = 1
    it = 0
    while it < nvary:
        it = it + 1
        ctl = r4(2.0) * auxiliary_routines.ran1() - r4(1.0)
        phil = auxiliary_routines.ran1() * r4(360.0)
        phi = auxiliary_routines.ran1() * r4(360.0)
        cz = r4(2.0) * auxiliary_routines.ran1() - r4(1.0)
        chi = auxiliary_routines.ran1() * r4(360.0)
        #            write(6,*) acos(ctl)*180.0/pi,phil,phi,acos(cz)*180.0/pi,chi
        for k in range(1, nexpan + 1):
            sumc = f32(0.0)
            sums = f32(0.0)
            spsumc = f32(0.0)
            spsums = f32(0.0)
            ll1 = int(l1ex[k])
            ll2 = int(l2ex[k])
            ll = int(lex[k])
            n1 = int(n1ex[k])
            n2 = int(n2ex[k])
            l1tmp = ll1
            l2tmp = ll2
            n1tmp = n1
            n2tmp = n2
            test = (l1only == 0 and l1tmp == 0) or l1only != 0
            test = test and ((l2only == 0 and l2tmp == 0) or l2only != 0)
            test = test and ((lonly == 0 and ll == 0) or lonly != 0)
            test = test and ((n1only == 0 and n1tmp == 0) or n1only != 0)
            test = test and ((n2only == 0 and n2tmp == 0) or n2only != 0)
            if test:
                phs = f32(1.0)
                phsn = f32(_m1pow(l1tmp + l2tmp + ll + n1tmp + n2tmp))

                # Normal algorithm for molecule 1 at the origin

                if identsign >= 0:
                    if m2fix > 0:
                        m2f = -l2tmp
                        m2l = l2tmp
                    else:
                        m2f = 0
                        m2l = 0
                    sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl,
                                                              phi, cz, chi, phs, sumc, sums)
                    # now reverse signs of n1,n2 if either is non-zero
                    if (n2 + n1) > 0:
                        n1tmp = -n1
                        n2tmp = -n2
                        phs = phs * phsn
                        sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil,
                                                                  ctl, phi, cz, chi, phs, sumc, sums)
                    if l2tmp == 0:  # Store the spatial terms separately
                        spsumc = sumc
                        spsums = sums
                else:
                    if m2fix > 0:
                        m2f = -l1tmp
                        m2l = l1tmp
                    else:
                        m2f = 0
                        m2l = 0

                    # Phase factor to reverse parity of centre to centre vector

                    phs = f32(_m1pow(ll))
                    sumc, sums = plotsharm_routines.sumoverm2(l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil, ctl,
                                                              phi, cz, chi, phs, sumc, sums)
                    # now reverse signs of n1,n2 if either is non-zero
                    if (n2 + n1) > 0:
                        n1tmp = -n1
                        n2tmp = -n2
                        phs = phsn * phs
                        sumc, sums = plotsharm_routines.sumoverm2(l2tmp, l1tmp, ll, n2tmp, m2f, m2l, n1tmp, phil,
                                                                  ctl, phi, cz, chi, phs, sumc, sums)
                    if l1tmp == 0:  # Store the spatial terms separately since reversing l1 and l2 will not lead to a
                        # spatial term.
                        spsumc = sumc
                        spsums = sums
            # now repeat for interchange of l1,n1,with l2,n2 for identical molecules
            # bearing in mind the phase convention on the h's - the case when l1=l2 AND
            # n1=n2 doesn't count here.
            #
            if ident == 0 and (ll1 != ll2):  # l1 will always be .ge. l2 (from the way they are calculated).
                phs = f32(_m1pow(ll1 + ll2))
                l1tmp = ll2
                l2tmp = ll1
                n1tmp = n2
                n2tmp = n1
                test = (l1only == 0 and l1tmp == 0) or l1only != 0
                test = test and ((l2only == 0 and l2tmp == 0) or l2only != 0)
                test = test and ((lonly == 0 and ll == 0) or lonly != 0)
                test = test and ((n1only == 0 and n1tmp == 0) or n1only != 0)
                test = test and ((n2only == 0 and n2tmp == 0) or n2only != 0)
                if test:
                    phsn = f32(_m1pow(l1tmp + l2tmp + ll + n1tmp + n2tmp))
                    if m2fix > 0:
                        m2f = -l2tmp
                        m2l = l2tmp
                    else:
                        m2f = 0
                        m2l = 0
                    sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil, ctl,
                                                              phi, cz, chi, phs, sumc, sums)
                    # now reverse signs of n1,n2 if either is non-zero
                    if (n2 + n1) > 0:
                        n1tmp = -n2
                        n2tmp = -n1
                        phs = phs * phsn
                        sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil,
                                                                  ctl, phi, cz, chi, phs, sumc, sums)
            zcoeffc[k] = sumc
            zcoeffspc[k] = spsumc
            if forcemirrorsym:
                zcoeffs[k] = 0.0
                zcoeffsps[k] = 0.0
            else:
                zcoeffs[k] = sums
                zcoeffsps[k] = spsums
            #               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
        for i in range(1, nr + 1):
            tsumrd = f32(0.0)
            tsumsp = f32(0.0)
            tsumto = f32(0.0)
            for k in range(1, nexpan + 1):
                term = datac[i, k] * zcoeffc[k] + datas[i, k] * zcoeffs[k]
                tsumto = tsumto + term
                if l2ex[k] == 0:
                    term = datac[i, k] * zcoeffspc[k] + datas[i, k] * zcoeffsps[k]
                    tsumsp = tsumsp + term
                    if l1ex[k] == 0:
                        tsumrd = tsumrd + term
            temprd[i] = tsumrd
            tempsp[i] = tsumsp
            tempto[i] = tsumto
        sumrd = float(kineticterm)
        sumsp = float(kineticterm)
        sumto = float(kineticterm)
        for i in range(1, nr + 1):
            factr = fact * rd[i] * rd[i]
            temp = temprd[i]
            term = r4(1.0) + temp
            if term > maxrd:
                maxrd = term
                maxrdr = rd[i]
            if term < minrd:
                minrd = term
                minrdr = rd[i]
            if term > 0.0:
                contribution = factr * (temp - term * logf(term))
            else:
                contribution = -factr
            sumrd = sumrd + float(contribution)
            grw8[i, 1] = grw8[i, 1] + sumrd
            temp = tempsp[i]
            term = r4(1.0) + temp
            if term > maxsp:
                maxsp = term
                maxspr = rd[i]
                thetalspmax = acosf(ctl) * r4(180.0) / pi
                philspmax = phil
                phimspmax = phi
                thetamspmax = acosf(cz) * r4(180.0) / pi
                chimspmax = chi
            if term < minsp:
                minsp = term
                minspr = rd[i]
                thetalspmin = acosf(ctl) * r4(180.0) / pi
                philspmin = phil
                phimspmin = phi
                thetamspmin = acosf(cz) * r4(180.0) / pi
                chimspmin = chi
            #                  if(term.lt.0.0) then
            #                     spnegative(i)=spnegative(i)-term
            #                  end if
            # Make correction for negative regions of g(r)
            if temp > 0.0 and spnegative[i] > 0.0:
                wtnegative = temp * spnegative[i] / term
                temp = temp + wtnegative
                spnegative[i] = spnegative[i] - wtnegative
                term = r4(1.0) + temp
            if term > 0.0:
                contribution = factr * (temp - term * logf(term))
            else:
                contribution = -factr
            sumsp = sumsp + float(contribution)
            grw8[i, 2] = grw8[i, 2] + sumsp
            temp = tempto[i]
            term = r4(1.0) + temp
            if term > maxto:
                maxto = term
                maxtor = rd[i]
                thetaltomax = acosf(ctl) * r4(180.0) / pi
                philtomax = phil
                phimtomax = phi
                thetamtomax = acosf(cz) * r4(180.0) / pi
                chimtomax = chi
            if term < minto:
                minto = term
                mintor = rd[i]
                thetaltomin = acosf(ctl) * r4(180.0) / pi
                philtomin = phil
                phimtomin = phi
                thetamtomin = acosf(cz) * r4(180.0) / pi
                chimtomin = chi
            #                  if(term.lt.0.0) then
            #                     tonegative(i)=tonegative(i)-term
            #                  end if
            # Make correction for negative regions of g(r)
            if temp > 0.0 and tonegative[i] > 0.0:
                wtnegative = temp * tonegative[i] / term
                temp = temp + wtnegative
                tonegative[i] = tonegative[i] - wtnegative
                term = r4(1.0) + temp
            if term > 0.0:
                contribution = factr * (temp - term * logf(term))
            else:
                contribution = -factr
            sumto = sumto + float(contribution)
            grw8[i, 3] = grw8[i, 3] + sumto
        entropyrd = entropyrd + sumrd
        entropysp = entropysp + sumsp
        entropyto = entropyto + sumto
        if it == nfreq * idiv(it, nfreq) or it == nvary:
            units.write(6, fmt("(5(e13.6,1x))", ctl, phil, phi, cz, chi))
            units.write(6, fmt("(a,i10,a,3(1x,e13.6))", "mapgrentropy> Entropy after ", it, " trials:",
                               entropyrd / float(f32(it)), entropysp / float(f32(it)), entropyto / float(f32(it))))
            units.write(10, fmt("(i10,6(1x,e13.6))", it, sumrd, sumsp, sumto, entropyrd / float(f32(it)),
                                entropysp / float(f32(it)), entropyto / float(f32(it))))
    units.write(10, lw())
    units.write(10, fmt('("Minimum values, radial, spatial total: ",3(1x,e12.5))', minrd, minsp, minto))
    units.write(10, fmt('("Radius values of minima, radial, spatial total: ",3(1x,e12.5))', minrdr, minspr, mintor))
    units.write(10, lw())
    units.write(10, fmt('("thetal, phil, phim, thetam, chim of spatial minimum: ",5(1x,e12.5))', thetalspmin,
                        philspmin, phimspmin, thetamspmin, chimspmin))
    units.write(10, fmt('("thetal, phil, phim, thetam, chim of total minimum: ",5(1x,e12.5))', thetaltomin,
                        philtomin, phimtomin, thetamtomin, chimtomin))
    units.write(10, lw())
    units.write(10, fmt('("Maximum values, radial, spatial total: ",3(1x,e12.5))', maxrd, maxsp, maxto))
    units.write(10, fmt('("Radius values of maxima, radial, spatial total: ",3(1x,e12.5))', maxrdr, maxspr, maxtor))
    units.write(10, lw())
    units.write(10, fmt('("thetal, phil, phim, thetam, chim of spatial maximum: ",5(1x,e12.5))', thetalspmax,
                        philspmax, phimspmax, thetamspmax, chimspmax))
    units.write(10, fmt('("thetal, phil, phim, thetam, chim of total maximum: ",5(1x,e12.5))', thetaltomax,
                        philtomax, phimtomax, thetamtomax, chimtomax))
    units.close(10)
    #
    # write out the results
    #
    units.open(10, "negativegr.txt", status="unknown")
    for i in range(1, nr + 1):
        units.write(10, fmt("(3(1x,e12.5))", rd[i], spnegative[i], tonegative[i]))
    units.close(10)
    for i in range(1, nr + 1):
        for it in range(1, 4):
            grw[i, it] = f32(grw8[i, it] / float(f32(nvary)))
            grwerr[i, it] = f32(grwerr8[i, it] / float(f32(nvary)))
        grw[i, 4] = grw[i, 2] - grw[i, 1]
        grwerr[i, 4] = 0.0
        grw[i, 5] = grw[i, 3] - grw[i, 2]
        grwerr[i, 5] = 0.0
    units.write(6, lw(grw[nr, 1], entropyrd / float(f32(nvary))))

    # set the axis title

    titler = "r[A]"
    titledat[1] = fstr(fmt("(10x,a,11x)", "S_radial"), 29)
    titledat[2] = fstr(fmt("(10x,a,10x)", "S_spatial"), 29)
    titledat[3] = fstr(fmt("(11x,a,12x)", "S_total"), 29)
    titledat[4] = fstr(fmt("(7x,a,7x)", "S_extra_spatial"), 29)
    titledat[5] = fstr(fmt("(7x,a,8x)", "S_extra_orient"), 29)
    ilen = auxiliary_routines.indexlastdot(symbols.fnameinp)
    ext = "s"
    fname = fstr(substr(symbols.fnameinp, 1, ilen) + ext, 256)
    fname = auxiliary_routines.w_int(nr, fname, 5, rd, grw, grwerr, titler, titledat, nxdim, 5)
    units.open(10, "ran1.dat", status="unknown")
    auxiliary_routines.saveran1(10)
    units.close(10)
