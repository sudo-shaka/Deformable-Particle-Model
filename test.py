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
C = Cell(1,0,0,1.17,50,1.0,0.1,1.0,0.1,0.1,0.2,1.0,0.0);

for i in range(1000000):
  C.UpdateDirectorDiffusion(0.005)
  C.UpdateEuler(0.005)
  if i % 100000 == 0:
    print(C.GetPerim(), C.GetArea())
    plt.scatter(C.X,C.Y);

plt.axis('equal')
plt.show()
