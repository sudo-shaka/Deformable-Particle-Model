import DPM
from matplotlib import pyplot as plt
from progressbar import progressbar
import imageio

def Euler():
  C = DPM.Cell(1,0,0,1.4,32,1.0,0.5,0.1,0.05,0.2,0.1,20.0,0)
#  C.l1 = [0.1]*32;
#  C.l2 = [0.05]*32;
#  C1.l1 = [0.1]*28;
#  C1.l2 = [0.05]*28;


  ncells = 7;
  mono = DPM.monolayer([C]*ncells,0.85);
  mono.disperse();
  count = 0;
  mono.Kc = 0.5;

  for i in progressbar(range(int(20000))):
    for ci in range(ncells):
      mono.Cells[ci].UpdateDirectorDiffusion(0.001)
    mono.UpdateEuler(0.001);
    if i % 1000 == 0:
      for ci in range(ncells):
        plt.scatter(mono.Cells[ci].X,mono.Cells[ci].Y)
      plt.axis('equal');
      plt.savefig('/tmp/'+str(count)+'.png')
      plt.clf();
      count +=1;

  with imageio.get_writer('/tmp/out.gif',mode='I') as writer:
    for i in range(count):
      filename = '/tmp/'+str(i)+'.png'
      image = imageio.imread(filename)
      writer.append_data(image)
  return mono


if __name__ == "__main__":
  Euler();
