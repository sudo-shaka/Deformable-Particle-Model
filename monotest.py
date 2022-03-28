try:
  import DPM
  from matplotlib import pyplot as plt
  from progressbar import progressbar
  import imageio
  from Plot import PlotDPM
except:
  print("You do not have the required dependencies to run this simulation");

def Euler():
  #initialize a cell with calA of 1.4, 32 vertecies, spring and director constants plus an area of 20.
  C = DPM.Cell(1,0,0,1.4,32,1.0,0.5,0.1,0.05,0.2,0.1,20.0,0)

  ncells = 7;
  #initialize a monolayer if n number cells and a packing fraction of 
  mono = DPM.monolayer([C]*ncells,0.85);
  #Disperse the monolayer using simple hard sphere repulsions
  mono.disperse();
  count = 0;
  #Set the spring constant for cell-cell interactions
  mono.Kc = 0.5;
  
  #begining time steps
  for i in progressbar(range(int(20000))):
    #for every cell change the direction its going
    for ci in range(ncells):
      mono.Cells[ci].UpdateDirectorDiffusion(0.001)
    #update the positions with a dt of 0.001
    mono.UpdateEuler(0.001);
    #plot every 1000 time steps
    if i % 1000 == 0:
      PlotDPM(mono)
      plt.axis('equal')
      plt.savefig('/tmp/'+str(count)+'.png')
      plt.clf();
      count +=1;
  
  #this is to simply save the image
  with imageio.get_writer('/tmp/out.gif',mode='I') as writer:
    for i in range(count):
      filename = '/tmp/'+str(i)+'.png'
      image = imageio.imread(filename)
      writer.append_data(image)
  #return the monolayer
  print('Image is saved to "/tmp/out.gif"')
  return mono


if __name__ == "__main__":
  Euler();
