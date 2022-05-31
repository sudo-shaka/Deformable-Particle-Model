try:
  from DPM import Cell
except:
  print('Please compile the module using "python3 setup.py build_ext --inplace"')
  exit(1)

try:
  from matplotlib import pyplot as plt
except:
  print('This demo required matplotlib pyplot to be installed. \nPlease install with "pip install matplotlib"')
  exit(1)

import imageio
print("This is an example of a single deformable particle crawling around")
#input order: index, centerX1,centerY1, CalA0, vertexNumber, Kl, Kb, Ka, velocity0, Dr1, Ds2, area0, director (psi)
C = Cell(0,0,2.0,75,1.0,0.05,1.0,0.0,0.0,0.0,1.0,0.0);

with imageio.get_writer('out.gif',mode='I') as w:
  for i in range(25000):
    C.UpdateDirectorDiffusion(0.001)
    C.UpdateEuler(1,0.0001)
    if i % 500 == 0:
      F = [abs(C.Fx[i])+abs(C.Fy[i]) for i in range(C.NV)]
      plt.scatter(C.X,C.Y,c=F,cmap='coolwarm');
      plt.axis('equal')
      plt.savefig('/tmp/'+str(i)+'.png')
      w.append_data(imageio.imread('/tmp/'+str(i)+'.png'))
      plt.clf()
