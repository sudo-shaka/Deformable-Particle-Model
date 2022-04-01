import DPM
from matplotlib import pyplot as plt
import numpy as np

def PlotDPM(monolayer):
    L = monolayer.BoxLength;
    plt.figure(figsize=(8,8))
    for ci in range(monolayer.NCELLS):
        X = []; Y = [];
        for vi in range(len(monolayer.Cells[ci].X)):
            x = monolayer.Cells[ci].X[vi]
            y = monolayer.Cells[ci].Y[vi]
            if monolayer.Cells[ci].X[vi] < 0:
                x = monolayer.Cells[ci].X[vi] + L*(int(abs(monolayer.Cells[ci].X[vi])/L)+1)
            if monolayer.Cells[ci].X[vi] > L:
                x = monolayer.Cells[ci].X[vi] - L*(int(abs(monolayer.Cells[ci].X[vi]-L)/L)+1)
            if monolayer.Cells[ci].Y[vi] < 0:
                y = monolayer.Cells[ci].Y[vi] + L*(int(abs(monolayer.Cells[ci].Y[vi])/L)+1)
            if monolayer.Cells[ci].Y[vi] > L:
                y = monolayer.Cells[ci].Y[vi] - L*(int(abs(monolayer.Cells[ci].Y[vi]-L)/L)+1)
            X.append(x)
            Y.append(y)
        F = [abs(monolayer.Cells[ci].Fx[i]) + abs(monolayer.Cells[ci].Fy[i]) \
            for i in range(len(monolayer.Cells[ci].X))] 
        plt.scatter(X,Y,c=F,cmap='coolwarm');
        #cx = np.mean(x)
        #cy = np.mean(y)
        #plt.scatter(cx,cy,color='black')
    plt.axis('equal')
    #plt.show()
