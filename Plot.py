import DPM
from matplotlib import pyplot as plt

#Simple script to rebox the periodic boundary conditions and then plot using matplotlib
def PlotDPM(monolayer):
    L = monolayer.BoxLength;
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
        plt.scatter(X,Y);

    plt.axis('equal')
    #plt.show()