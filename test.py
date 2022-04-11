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


print("This is an example of a single deformable particle crawling around")
#input order: index, centerX1,centerY1, CalA0, vertexNumber, Kl, Kb, Ka, velocity0, Dr1, Ds2, area0, director (psi)
C = Cell(0,0,1.17,40,1.0,0.05,1.0,0.1,0.05,0.05,1.0,0.0);

for i in range(1000000):
  C.UpdateDirectorDiffusion(0.001)
  C.UpdateEuler(1,0.001)
  if i % 200000 == 0:
    plt.scatter(C.X,C.Y);

plt.axis('equal')
plt.show()
