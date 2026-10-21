c     Driver for tests/library/test_setup_symbols_generic.py: calls every routine of
c     fortran_src/setup_symbols_generic.f and dumps the symbol table after each call.
c     Build and run (from tests/library/data/setup_symbols, where the .symbols files are):
c
c       F=../../../../fortran_src
c       gfortran -O2 -ffixed-line-length-none -fno-backslash -w -I$F -o /tmp/ssg
c    $     ../setup_symbols_driver.f $F/setup_symbols_generic.f $F/parse.f
c       /tmp/ssg > ../setup_symbols_reference.txt
c
c     The Python side (test_setup_symbols_generic.py) writes the same dump.
      program ssgdriver
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      include 'commands.inc'
      character*256 auxname
      character*20 cases(16)
      integer ierr, k
      data cases /'good extra words', 'nondist', 'noeol', 'noprompt',
     $     'badprompt', 'noext', 'nodataprompt', 'nodptext',
     $     'blankinside', 'badkey', 'eofinside', 'nodefault', 'notype',
     $     'notext', 'missing', '  good'/

      call hdr('changeato')
      call setup_symbols_changeato()
      call dump()
      call hdr('partials')
      call setup_symbols_partials()
      call dump()
      call hdr('coord')
      call setup_symbols_coord()
      call dump()
      call hdr('clusters')
      call setup_symbols_clusters()
      call dump()
      call hdr('rings')
      call setup_symbols_rings()
      call dump()
      call hdr('chains')
      call setup_symbols_chains()
      call dump()
      call hdr('lsi')
      call setup_symbols_lsi()
      call dump()
      call hdr('sharm')
      call setup_symbols_sharm()
      call dump()
      call hdr('sdf')
      call setup_symbols_sdf()
      call dump()
      call hdr('sdfcube')
      call setup_symbols_sdfcube()
      call dump()
      call hdr('mapgrentropy')
      call setup_symbols_mapgrentropy()
      call dump()
      call hdr('plot2d')
      call setup_symbols_plot2d()
      call dump()
      call hdr('splot2d')
      call setup_symbols_splot2d()
      call dump()
      call hdr('plot3d')
      call setup_symbols_plot3d()
      call dump()
      call hdr('plot3djmol')
      call setup_symbols_plot3djmol()
      call dump()
      call hdr('mapgr')
      call setup_symbols_mapgr()
      call dump()
      call hdr('triangles')
      call setup_symbols_triangles()
      call dump()
      call hdr('torangles')
      call setup_symbols_torangles()
      call dump()
      call hdr('voids')
      call setup_symbols_voids()
      call dump()
      call hdr('fluctuations')
      call setup_symbols_fluctuations()
      call dump()
      call hdr('nwts')
      call setup_symbols_nwts()
      call dump()
      call hdr('xwts')
      call setup_symbols_xwts()
      call dump()
      call hdr('bragginten')
      call setup_symbols_bragginten()
      call dump()
c     sharm again: its read of symbol_value(nsymbol+1) now sees bragginten's leftovers
      call hdr('sharm')
      call setup_symbols_sharm()
      call dump()
      do k=1,16
         call hdr('generic '//trim(cases(k)))
         auxname=cases(k)
         ierr=99
         call setup_symbols_generic(auxname,ierr)
         write(6,'(a,i4)') 'ssg_ierr ',ierr
         write(6,'(a,3i6)') 'commands ',nwords,ncf(1),ncl(1)
         call dump()
      end do
      end

      subroutine hdr(name)
      character*(*) name
      write(6,'(a,a)') '=== ',name
      end

      subroutine dump()
      include 'epsr_dimension.inc'
      include 'symbols.inc'
      integer i
      write(6,'(a,6i6)') 'n ',nsymbol,nsymboldata,nsymboluptodata,
     $     nsymbolplus,nfiles,irefrho
      write(6,'(a,i6,a)') 'prompt ',ilenprompt,' |'//trim(prompt)//'|'
      write(6,'(a,i6,a)') 'dataprompt ',ilendataprompt,
     $     ' |'//trim(dataprompt)//'|'
      write(6,'(a)') 'dataprompttext |'//trim(dataprompttext)//'|'
      write(6,'(a,i6,a)') 'extinp ',ilenextinp,' |'//trim(extinp)//'|'
      do i=1,nsymbol
         write(6,'(i4,1x,i3,1x,a,a,a,a,a,a,a,a,a,a)') i,iflag(i),
     $        symbol(i),'|',symbol_type(i),'|',trim(symbol_value(i)),
     $        '|',trim(symbol_default(i)),'|',trim(symbol_text(i)),'|'
      end do
      end
