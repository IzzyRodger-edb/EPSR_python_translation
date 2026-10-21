"""Translation of fortran_src/bondnumbers.f."""

from ..common import clustercom, molecom, molecom_mod
from ..trace import routine


@routine
def bondnumbers():
    """bondnumbers.f:bondnumbers -- renumber bond atoms relative to each molecule's first atom.

    to determine the correct atom numbers to assign to bonds in the input .ato
    file. This is to correct for the possibility that the atom numbers used to
    describe the intramolecular bonds may have got out sync. with the atom
    numbering in this file

    The method is to step through all the atom numbers used to define bonds.
    The lowest of these is assumed to refer to the first atom in the molecule.
    The bond atom numbers are then revised to be in the correct relationship
    with the actual atom number of the first atom in the molecule
    """
    # use molecom_mod
    naf = clustercom.naf
    nal = clustercom.nal
    nbond = molecom.nbond
    iabond = molecom.iabond
    imageatom = molecom.imageatom
    ig = 0
    for im in range(1, clustercom.nmol + 1):
        namin = 0
        for ia in range(naf[im], nal[im] + 1):
            if nbond[ia] > 0:
                for ib in range(1, nbond[ia] + 1):
                    nt = int(iabond[ia, ib])
                    if namin == 0:
                        namin = nt
                    if nt < namin:
                        namin = nt
        #     write(6,*) im,itype(naf(im)),naf(im),namin
        #
        # reset the bond atom numbers if necessary
        #
        nafim = int(naf[im])
        if namin != 0 and namin != nafim:
            for ia in range(naf[im], nal[im] + 1):
                if imageatom[ia] > 0:
                    imageatom[ia] = int(imageatom[ia]) - namin + nafim
                if nbond[ia] > 0:
                    for ib in range(1, nbond[ia] + 1):
                        nt = int(iabond[ia, ib]) - namin
                        iabond[ia, ib] = nafim + nt
        if clustercom.ngrot[im] > 0:
            for ih in range(1, clustercom.ngrot[im] + 1):
                ig = ig + 1
                #
                # reset the internal rotation atom numbers if necessary
                #
                if namin != 0 and namin != nafim:
                    labaxi = molecom_mod.labaxi
                    labrot = molecom_mod.labrot
                    labaxi[ig, 1] = int(labaxi[ig, 1]) - namin + nafim
                    labaxi[ig, 2] = int(labaxi[ig, 2]) - namin + nafim
                    for j in range(molecom_mod.rotf[ig], molecom_mod.rotl[ig] + 1):
                        labrot[j] = int(labrot[j]) - namin + nafim
                    if molecom_mod.rotgrp[ig] == "DIH":
                        for k in range(1, molecom.ndihedralang[ig] + 1):
                            for j in range(1, 3):
                                molecom.labrotdih[j, k, ig] = int(molecom.labrotdih[j, k, ig]) - namin + nafim
    return
