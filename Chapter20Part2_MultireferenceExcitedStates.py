# This Python script uses high-level multireference perturbation theory
# calculations to compute the first 12 excited states of beryllium atom.
# Calculations use a four-orbital, two-electron active space to treat the two
# valence electrons, combined with N-electron valence perturbation theory to
# treat remaining dynamical correlation. 
#
from pyscf import gto,scf,mcscf,mrpt
import numpy 

m=gto.Mole(atom='Be',basis='aug-cc-pvqz',charge=0,spin=0) 
m.build()
mf=scf.UHF(m)
mf.kernel()
wts=numpy.array((.3, .1,.1,.1, .1,.1,.1, .1, .1,.1,.1,  .1))
wts=wts/numpy.sum(wts)

mc=mcscf.CASSCF(mf,4,2).state_average_(weights=numpy.asarray(wts))
mc.kernel()
orbs=mc.mo_coeff


mc2=mcscf.CASCI(mf,4,2)
mc2.fcisolver.nroots=12
mc2.kernel(orbs) 
ncas = mc2.ncas
ncore = mc2.ncore
actsl = slice(ncore,ncore+ncas)
actcsl = slice(0,ncore+ncas)
csl = slice(0,ncore)

for i in range(mc2.fcisolver.nroots):
  print('PROCESSING STATE ',i)
  e_corr = mrpt.NEVPT(mc2,root=i).kernel()
  e_tot = mc2.e_tot[i] + e_corr
  print('E',mc2.e_tot[i],e_tot)


