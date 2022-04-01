def OutputMonolayerToFile(monolayer,outfilename):
  PHI = monolayer.GetPackingFraction()
  NCELLS = len(monolayer.Cells)
  L = monolayer.BoxLength
  with open(outfilename,'a') as out:
    out.write("START\nPHI,"+str(PHI))
    out.write(",NCELLS,"+str(NCELLS))
    out.write(",L,"+str(L)+'\n')
    for ci in range(NCELLS):
      out.write("CELLINFO,"+str(ci)+ ",NV,"+str(len(monolayer.Cells[ci].X)+1)\
          +",CalA,"+str(monolayer.Cells[ci].GetShapeParameter())+'\n')
      out.write("X,"+str(monolayer.Cells[ci].X)[1:-1]+"\n")
      out.write("Y,"+str(monolayer.Cells[ci].Y)[1:-1]+"\n")
      out.write("FX,"+str(monolayer.Cells[ci].Fx)[1:-1]+"\n")
      out.write("FY,"+str(monolayer.Cells[ci].Fy)[1:-1]+"\n")

    out.close()

def OutputCellToFile(Cell,outfilename):
  with open(outfilename,'a') as out:
    out.write("Cell: ,NV,"+str(len(Cell.X+1)))
    out.write("X,"+str(Cell.X)[1:-1]+'\n')
    out.write("Y,"+str(Cell.Y)[1:-1]+'\n')
    out.write("FX,"+str(Cell.Fx)[1:-1]+'\n')
    out.write("FX,"+str(Cell.Fy)[1:-1]+'\n')
    out.close()

def ReadCellFromFile(infilename):
  try:
    import DPM
  except:
    print("Error! DPM module not installed")
  with open(infilename,'r') as f:
    lines = f.readlines()
    count = 0
    CellArr = []
    for l in lines:
      data = l.split(",")
      if l[:4] == "Cell":
        if count > 0:
          CellArr.append(C)
        nv = int(data[3])
        C = DPM.Cell(nv)
        count += 1
      elif l[:1] == "X":
        C.X = [float(i) for i in data[1:]]
      elif l[:1] == "Y":
        C.Y = [float(i) for i  in data[1:]]
      elif l[:2] == "FX":
        C.Fx = [float(i) for i in data[1:]]
      elif l[:2] == "FY":
        C.Fx = [float(i) for i in data[1:]]

    f.close()
    return CellArr



def ReadMonolayerFromFile(infilename):
  try:
    import DPM
  except:
    print("Error! DPM module not installed")

  countm = 0
  monolist = []; PHIlist = []; NVlist = [];
  with open(infilename,'r') as f:
    lines = f.readlines()
    for l in lines:
      data = l.split(',')
      if l[:5] == "START":
        if countm > 0:
          mono = DPM.monolayer(CellArr,PHIlist[countm-1]);
          mono.BoxLength = L
          monolist.append(mono)
        countm += 1
        countc = 0
        CellArr = [];
      elif l[:3] == "PHI":
        PHIlist.append(float(data[1]))
        L = float(data[5])
      elif l[:8] == "CELLINFO":
        if countc > 0:
          CellArr.append(C)
        nv = int(data[3])
        C = DPM.Cell(nv)
        countc += 1
      elif l[:1] == "X":
        x = [float(i) for i in data[1:]]
        C.X = x
      elif l[:1] == "Y":
        y = [float(i) for i in data[1:]]
        C.Y = y
      elif l[:2] == "FX":
        Fx = [float(i) for i in data[1:]]
        C.Fx = Fx
      elif l[:2] == "FY":
        Fy = [float(i) for i in data[1:]]
        C.Fy = Fy

  f.close()
  return monolist
