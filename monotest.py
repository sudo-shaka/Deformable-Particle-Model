#This code runs the example seen in the README
from platform import platform


try:
  import DPM
  from matplotlib import pyplot as plt
  from progressbar import progressbar
  import imageio
  from Plot import PlotDPM
  import DPMIO
except:
  print("You do not have the required dependencies to run this simulation");

def Euler():

  mono = DPMIO.ReadMonolayerFromParams('input.csv')
  #begining time steps
  dt = 0.001; nsteps = 100; nout = 50;
  for i in progressbar(range(nout)):
    PlotDPM(mono)
    plt.axis('equal')
    plt.savefig('/tmp/'+str(i)+'.png')
    plt.close();
    mono.UpdateEuler(nsteps,dt);

  #this is to simply save the image
  with imageio.get_writer('/tmp/out.gif',mode='I') as writer:
    for i in range(nout):
      filename = '/tmp/'+str(i)+'.png'
      image = imageio.imread(filename)
      writer.append_data(image)
  #return the monolayer
  print('Image is saved to "/tmp/out.gif"')
  return mono

if __name__ == "__main__":
  Euler();
