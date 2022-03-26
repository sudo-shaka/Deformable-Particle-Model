from DPM import Cell
from matplotlib import pyplot as plt

C = Cell(1,0,0,1.17,50,1.0,0.1,1.0,0.1,0.1,0.2,1.0,0.0);

for i in range(1000000):
  C.UpdateDirectorDiffusion(0.005)
  C.UpdateEuler(0.005)
  if i % 100000 == 0:
    print(C.GetPerim(), C.GetArea())
    plt.scatter(C.X,C.Y);

plt.axis('equal')
plt.show()
