# This script generates Hartree-Fock, PBE0 DFT, MP2, CCSD, CASSCF, and NEVPT2
# dissociation curves for spin-symmetry-restricted ground-state H2 molecule.
# The code uses the Hartree-Fock guess density matrix from short bond lengths
# to do the self-consistent field calculations for long bond lengths. Outputs
# are prefaced with "RES", and the "negative" bond length is used to denote
# energies of isolated H atom. 
# 
from pyscf import scf,gto,dft,mp,cc,ci, mcscf, mrpt 
import os 
import numpy ,scipy, sys , math 

xc='.25*hf+.75*pbe,pbe'
fb='aug-cc-pvqz'
nat=2
Pold=None

## Loop over molecules 
print("# Bond lengths in Angstrom, total energies in Hartree, negative bond length is isolated H atom \n")
print("# Bondlength HF MP2 CCSD CAS(2,2)SCF NEVPT2, DFT \n")
for r in (-1,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5.,6.):

  # Build the molecule
  geom=''
  for i in range(nat):
    geom = geom + 'H  %.2f 0.0 0.0 \n'%(i*r)
  charge=0
  spin=0
  if(r<0):
    charge=0
    spin=1
    geom='H'
    Pold=None

  m=gto.Mole(atom=geom,basis=fb,charge=charge,spin=spin) 
  m.build() 
  (Na,Nb)=m.nelec
  NAO=m.nao

  mf=scf.RHF(m) 
  if(Pold is not None):
    mf.kernel(dm0=Pold)
  else: 
    mf.kernel()
  EHF = mf.e_tot
  P = mf.make_rdm1() 
  mf2=mp.MP2(mf)
  mf2.run() 
  EMP2 = mf2.e_tot 
  mc = cc.CCSD(mf)
  mc.run()
  ECC = mc.e_tot
  ECAS = EHF
  ENEV = EHF 
  if(r>0):
    mcas = mcscf.CASSCF(mf,2,2)
    mcas.run()
    ECAS = mcas.e_tot
    EC=mrpt.NEVPT(mcas).kernel()
    ENEV = ECAS + EC 
  md=dft.RKS(m,xc=xc)
  md.kernel(dm0=P)
  EDFT=md.e_tot
  print('RES %.3f %.6f %.6f %.6f %.6f %.6f %.6f' % (r,EHF,EMP2,ECC,ECAS,ENEV,EDFT))
  
