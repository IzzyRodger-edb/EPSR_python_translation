"""Translation of fortran_src/mapgrquad.f: ``program standalone`` and
``subroutine run_mapgr`` -- this version of MAP3D puts out the g(r) in specified
directions, from the spherical-harmonic coefficients written by sharm/sdf.

``setup_symbols_mapgr`` is in setup_symbols_generic.f (library/setup_symbols_generic.py,
imported inside the routine: it is translated by another work package).  smoosharm and
sumoverm2 are in plotsharm_routines.f.

``external cosd,sind,...``: the cosd used here is EPSR's own (auxiliary_routines.cosd);
sumoverm2 (plotsharm_routines.f) uses gfortran's intrinsics.
"""

from ..common import commands, symbols
from ..fortran import PI4ATAN, Ch, FortranEOF, FortranIOError, FortranStop, f32, findex, fmt, fstr, int_array
from ..fortran import iostat_of, isign, len_trim, r4, read_internal, real_array, char_array, substr, units
from ..library import auxiliary_routines
from ..library import file_and_directory_commands
from ..library import lngamma_mod
from ..library import parse
from ..library import plotsharm_routines
from ..library import read_and_save_files
from ..trace import routine

# parameter(mexpan=3000,mval=10,mpixel=80,mpixel2=(2*mpixel+1),mpixel3=mpixel2**3)
mexpan = 3000
mval = 10
mpixel = 80
mpixel2 = 2 * mpixel + 1
mpixel3 = mpixel2 ** 3

_F92 = "(1x,10(1x,i5))"
_F359 = "(10x,'view_',i3.3)"


def main():
    """mapgrquad.f:standalone -- program to perform the EPSR setup command.

    2009/3/30 Attempt to make this a standalone program
    """
    # Get the command from the command line arguments
    file_and_directory_commands.get_command_from_args()
    # Run the command
    run_mapgr()
    raise FortranStop(0)


def _m1pow(n):
    """(-1)**n for an integer n (integer result; also right for negative n)."""
    return 1 if n % 2 == 0 else -1


def _read_index(temptitle, ifirst, lentemptitle, ch):
    """index(temptitle(ifirst:lentemptitle),ch)"""
    return findex(substr(temptitle, ifirst, lentemptitle), ch)


@routine
def run_mapgr():
    """mapgrquad.f:run_mapgr -- g(r) of the centres, and g(r,w) along the specified
    directions/relative orientations ("views"), written to <input>.g01."""
    from ..library import setup_symbols_generic

    # use lngamma_mod
    # dimension datac(500,mexpan),datas(500,mexpan),rd(500)
    # 1,lvl(3,mval),nl(3),ml(3),zcoeffc(mexpan),zcoeffs(mexpan),grw(500,mexpan),grwerr(500,mexpan)
    # 1,th(50),ph(50),xl(2),yl(2),temp(500)
    # 2,out(0:mpixel3),xx(201),yy(201),l1ex(mexpan)
    # 3,l2ex(mexpan),lex(mexpan),n1ex(mexpan),n2ex(mexpan)
    # *,rdcoor(50),thcoor(50),phcoor(50),radbal(50),rgbbal(3,50)
    # *,tempc(500),temps(500)
    # (lvl, nl, ml, th, ph, xl, yl, out, xx, yy, rdcoor, thcoor, phcoor, radbal, rgbbal are not used)
    datac = real_array(500, mexpan)
    datas = real_array(500, mexpan)
    rd = real_array(500)
    zcoeffc = real_array(mexpan)
    zcoeffs = real_array(mexpan)
    grw = real_array(500, mexpan)
    grwerr = real_array(500, mexpan)
    temp = real_array(500)
    l1ex = int_array(mexpan)
    l2ex = int_array(mexpan)
    lex = int_array(mexpan)
    n1ex = int_array(mexpan)
    n2ex = int_array(mexpan)
    tempc = real_array(500)
    temps = real_array(500)
    # character*256 infile,fname,fnames,fnameout
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
    cosd = auxiliary_routines.cosd
    # integer idot,nfreq,ident,indentsign,leninfile,lenfname
    # integer nsmoord,nsmoosp,nsmooor,grtype,nrminsmoo
    #
    # setup the mapgrentropy symbols
    #
    setup_symbols_generic.setup_symbols_mapgr()
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
    pi = PI4ATAN  # noqa: F841  4.0*atan(1.0)
    nxdim = 500
    nfdim = mexpan
    #
    # initialise log(factorial)
    #
    lngamma_mod.fctrln(250)
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
    iref = iref + 1
    nviews = 0
    nviews = read_internal(symbols.symbol_value[iref]).int(nviews)
    # Add an extra view - the centres distribution function
    ituse = 1
    for i in range(1, nr + 1):
        grw[i, ituse] = r4(1.0) + datac[i, 1]
    titledat[ituse] = fstr(fmt("(8x,a,9x)", "g_centres(r)"), 29)
    # NB: grtype, thetal, phil, phim, thetam, chim are not initialised in the Fortran; they are
    # only used after a successful read in the reference cases (0 as in the Rust mapgr.rs)
    grtype = 0
    thetal = f32(0.0)
    phil = f32(0.0)
    phim = f32(0.0)
    thetam = f32(0.0)
    chim = f32(0.0)
    # calculate g(r,w) at specified theta and phi values
    for it in range(1, nviews + 1):
        iref = iref + 1
        try:
            grtype = read_internal(symbols.symbol_value[iref]).int(grtype)
            itest = 0
        except (FortranEOF, FortranIOError) as e:
            itest = iostat_of(e)
        iref = iref + 1
        try:
            r = read_internal(symbols.symbol_value[iref])
            thetal = r.real(thetal)
            phil = r.real(phil)
            phim = r.real(phim)
            thetam = r.real(thetam)
            chim = r.real(chim)
            itest = 0
        except (FortranEOF, FortranIOError) as e:
            itest = iostat_of(e)
        l1only = 1
        l2only = 1
        lonly = 1
        n1only = 1
        n2only = 1
        m2fix = 1
        if grtype < 3:
            l2only = 0
            n2only = 0
            m2fix = 0
        if grtype < 2:
            l1only = 0
            lonly = 0
            n1only = 0
        if itest == 0:
            ctl = cosd(thetal)
            phi = phim
            cz = cosd(thetam)
            chi = chim
            #            write(6,*) acos(ctl)*180.0/pi,phil,phi,acos(cz)*180.0/pi,chi
            for k in range(1, nexpan + 1):
                sumc = f32(0.0)
                sums = f32(0.0)
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
                        # now reverse signs of n1,n2 if either are non-zero, not including when both are zero
                        if (n2 + n1) > 0:
                            n1tmp = -n2
                            n2tmp = -n1
                            phs = phs * phsn
                            sumc, sums = plotsharm_routines.sumoverm2(l1tmp, l2tmp, ll, n1tmp, m2f, m2l, n2tmp, phil,
                                                                      ctl, phi, cz, chi, phs, sumc, sums)
                zcoeffc[k] = sumc
                if forcemirrorsym:
                    zcoeffs[k] = 0.0
                else:
                    zcoeffs[k] = sums
                #               write(6,*) ll1,ll2,ll,n1,n2,phil,ctl,phi,cz,sumc,sums
            for i in range(1, nr + 1):
                sumt = f32(0)
                for k in range(1, nexpan + 1):
                    term = datac[i, k] * zcoeffc[k] + datas[i, k] * zcoeffs[k]
                    sumt = sumt + term
                temp[i] = sumt + r4(1.0)
            ituse = ituse + 1
            for i in range(1, nr + 1):
                grw[i, ituse] = temp[i]
                grwerr[i, ituse] = 0.0
            # Write the title for this distribution
            titledat[ituse] = fstr(fmt(_F359, ituse - 1), 29)
    #
    # write out the results
    #

    # set the axis title

    titler = "r[A]"
    # Generate name of output file
    idot = auxiliary_routines.indexlastdot(symbols.fnameinp)
    if idot > 0:
        fnameout = fstr(substr(symbols.fnameinp, 1, idot) + "g", 256)
    else:
        fnameout = fstr(substr(symbols.fnameinp, 1, len_trim(symbols.fnameinp)) + ".g", 256)
    units.write(6, fmt("(a,1x,a)", substr(symbols.fnameinp, 1, len_trim(symbols.fnameinp)),
                       substr(fnameout, 1, len_trim(fnameout))))
    # NB: grwerr(:,1) is never set in the Fortran (stack garbage in the error column of
    # g_centres); here it is 0.  compare_outputs.py masks those columns.
    fnameout = auxiliary_routines.w_int(nr, fnameout, ituse, rd, grw, grwerr, titler, titledat, nxdim, nfdim)
