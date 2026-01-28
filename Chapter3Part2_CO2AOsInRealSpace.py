# This Python script generates the values of each atomic orbital in STO-3G CO2,
# plotted along the bond axis (chosen as z). Note that the px and py atomic
# orbitals on each atom are zero at every point along the z axis. 
from pyscf import gto,scf
import numpy 
m=gto.Mole(atom='O 0.0 0.0 -1.16; C 0.0 0.0 0.0; O 0.0 0.0 1.16',basis='sto-3g')
m.build()
zs0=[]
for i in range(-200,200):
  zs0.append(0.02*i)
zs=numpy.array(zs0)
xs=numpy.zeros_like(zs)
coords=numpy.array((xs,xs,zs)).T
aos = m.eval_gto("GTOval_cart",coords)
#print('teh aos shape',aos.shape)
(nr,nao)=aos.shape
print("# z chi1(z) chi2(z).... \n")
for r in range(nr):
  st='%.3f '%(zs[r])
  for i in range(nao):
    st=st+'%.4f '%(aos[r,i])
  print(st)

