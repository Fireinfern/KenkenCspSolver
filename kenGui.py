from PyQt5.QtWidgets import *
from ortools.sat.python import cp_model
import itertools
import numpy as np
import inspect
import os

#Suma
#Pretty self explanatory
def addSumConstraint(x, n, celdas, model):
  model.Add(sum(celdas) == x)

## Multiplicacion utilizando una variable auxiliar
# x es nuestro resultado, n es la dimension del kenken, celdas es el arreglo de celdas que ocupa la operacion y model es el modelo de CSP
def addMultiplication(x, n, celdas, model):
  #model.AddDecisionStrategy(celdas, cp_model.CHOOSE_)
  nC = len(celdas)
  #Si solo son dos variables que se multiplican, no hay necesidad de crear un auxiliar
  #Por ello se usa el Multiplication ecuality con las dos celdas
  if nC == 2:
    model.AddMultiplicationEquality(x,celdas)
  else:
    #para el resto se crea un arreglo auxiliar que guarda el resultado de ambas
    aux =[]
    for i in range(nC):
      if i == 0:
        aux += [model.NewIntVar(1,1000, 'aux'+str(i))]
        model.AddMultiplicationEquality(aux[i],celdas[0:2])
      elif i > 1:
        aux += [model.NewIntVar(1,1000, 'aux'+str(i))]
        model.AddMultiplicationEquality(aux[i-1],[aux[i-2],celdas[i]])
    #al final se iguala el utimo axiliar a el resultado x
    model.Add(aux[len(aux)-1] == x)


def addResta(x,n,celdas, model, number):
  p = model.NewBoolVar('rest'+str(number))
  model.Add((celdas[0]-celdas[1] == x)).OnlyEnforceIf(p)
  model.Add((celdas[1]-celdas[0] == x)).OnlyEnforceIf(p.Not())
  #model.AddBoolXOr([(celdas[0]-celdas[1] == x),(celdas[1]-celdas[0] == x)])

def addDivision(x, n, celdas, model, number):
  q = model.NewBoolVar('div' + str(number))
  model.Add(celdas[0]*x==celdas[1]).OnlyEnforceIf(q)
  model.Add(celdas[1]*x==celdas[0]).OnlyEnforceIf(q.Not())
  #model.AddBoolXOr([b,c])

#Clase de Agrupación: guarda como estructura una operación de números con resultado
class Agrupation:
  def __init__(self,result,operation,positions):
      self.result = result        #Resultado de la Operación - int
      self.operation = operation      #Operación Aritmética - cadena string de len=1
      self.positions = positions[:]       #Array de Posiciones - [(x1,y2),...(Xn,Yn)]
      self.numbers = [0]*len(self.positions)      #Array de números soluciones

  def getPositions(self):
      return self.positions[:]

  def getNumbers(self):
      return self.numbers[:]
  
  def getResult(self):
      return self.result
  
  def getOperation(self):
      return self.operation


#Clase de apoyo para conseguir la posición de la agrupación
def returnPosition(listam,dimension,criteria):
    dum_lista = listam[:]
    result = None
    contador = 0
    for agrupacion in dum_lista:
      if dimension == 0:
        if agrupacion.getResult() == criteria: result = criteria
      elif dimension == 1:
        result = agrupacion.getOperation().index(criteria) if criteria in agrupacion.getOperation() else None
      elif dimension == 2:
        result = agrupacion.getPositions().index(criteria) if criteria in agrupacion.getPositions() else None
      if result != None: return contador
      contador = contador + 1
    return result

#Clase de la Lista de Agrupaciones
class AgrupationList:
  #Iniciar con una lista vacía o ya creada
  def __init__(self,agrupation_list = None, size = None):
    if(agrupation_list == None):
      self.agrupation_list = []
    else:
      self.agrupation_list = agrupation_list[:]
    if(size == None):
      self.size = 0
    else:
      self.size = size

  def addAgrupacion(self,agrupation):
    temp_agrupation = Agrupation(agrupation.result,agrupation.operation,agrupation.positions)
    self.agrupation_list.append(temp_agrupation)
  
  #Encontrar la agrupación según un criterio de búsqueda: Resultado, Operación, Coordenada de la celda
  def getAgrupacion(self,criteria):
    dimension = None
    try:
      dimension = len(criteria)
    except:
      dimension = 0
    pos_reference = returnPosition(self.agrupation_list[:], dimension, criteria)
    print(pos_reference)
    return self.agrupation_list[pos_reference]
  
  def getLength(self):
    return len(self.agrupation_list)
  
  def getAllAgrupations(self):
    return self.agrupation_list[:]
  
  def getSize(self):
    return self.size

def getMax(actual_max, val1, val2):
  if val1 > val2:
    if actual_max > val1:
       return actual_max
    else:
      return val1
  else:
    if actual_max > val2:
       return actual_max
    else:
      return val2

def FileReader(path_file):
  try:
    file = open(path_file,'r')
  except:
      return
  temp_list = []
  temp_max = -1
  for line in file:
    result_operation = line.split(',',2)
    pos = result_operation[2].split('\n',1)
    pos = pos[0].replace('[[', '')
    pos = pos.replace(']]','')
    pos = pos.replace('[','')
    pos = pos.replace(']','')
    pos = pos.split(',')
    array_pos = []

    for i in range(0,len(pos),2):
      temp = (int(pos[i]),int(pos[i+1]))
      array_pos.append(temp)
      temp_max = getMax(temp_max,temp[0],temp[1])

    temp_agrupation = Agrupation(int(result_operation[0]),result_operation[1],array_pos)

    temp_list.append(temp_agrupation)
  temp_max = temp_max + 1
  
  return (temp_list,temp_max)    
########


def CSP_kenkenSolver(data_path):
  temp = FileReader(data_path)
  agrupacion = AgrupationList(temp[0],temp[1])
  KenKen = agrupacion.getAllAgrupations()
  print(KenKen)
  model = cp_model.CpModel()
  ken = []
  n = agrupacion.getSize()
  for i in range(n):
    fila = []
    for j in range(n):
      fila += [model.NewIntVar(1, n, 'x'+str(i)+str(j))]
    ken += [fila]
  for i in range(len(KenKen)):
    celdas = []
    for j in range(len(KenKen[i].getPositions())):
      celdas += [ken[KenKen[i].getPositions()[j][0]][KenKen[i].getPositions()[j][1]]]
    if(KenKen[i].getOperation() == '+'):
      model.AddDecisionStrategy(celdas, cp_model.CHOOSE_FIRST,cp_model.SELECT_MIN_VALUE)
      addSumConstraint(KenKen[i].getResult(),n,celdas,model)
    if(KenKen[i].getOperation() == '-'):
      addResta(KenKen[i].getResult(),n,celdas,model,i)
    if(KenKen[i].getOperation() == '*'):
      addMultiplication(KenKen[i].getResult(),n,celdas,model)
    if(KenKen[i].getOperation() == '/'):
      addDivision(KenKen[i].getResult(),n,celdas,model,i)
    if(KenKen[i].getOperation() == ''):
      model.Add(celdas[0] == KenKen[i].getResult())
      if(KenKen[i].getResult() > n//2):
        model.AddDecisionStrategy(celdas,cp_model.CHOOSE_FIRST,cp_model.SELECT_UPPER_HALF)
      else:
        model.AddDecisionStrategy(celdas,cp_model.CHOOSE_FIRST,cp_model.SELECT_LOWER_HALF)
  #Los Numeros en las filas son diferentes
  for i in range(n):
    model.AddAllDifferent(ken[i])
  #los numeros de las columnas deben ser diferentes
  for i in range(n):
    a = []
    for j in range(n):
        a += [ken[j][i]]
    model.AddAllDifferent(a)
  solver = cp_model.CpSolver()
  status = solver.Solve(model)
  final =[]
  if status == cp_model.FEASIBLE:
    for i in range(n):
      fila = []
      for j in range(n):
          fila +=[solver.Value(ken[i][j])]
      final+=[fila]
  return final

####################################################################################################
####################################################################################################
##App Creation######################################################################################
####################################################################################################
####################################################################################################
app = QApplication([])
MainTab = QTabWidget()
MainTab.setWindowTitle('Kenken Csp Solver')
CreateWindow = QWidget()
mainLayout = QVBoxLayout()
layout = QHBoxLayout()
Label_Ken = QLabel('Write the kenken size to be solve:')
nCells = QSpinBox()
nCells.setMaximum(10)
nCells.setMinimum(1)
listaRes = QListWidget()
listaOp = QListWidget()
def CreateGrid():
    nCells.setReadOnly(True)
    CreateButton.setDisabled(True)
    MainTab.setGeometry(30,30,800,600)
    GridLayout = QHBoxLayout()
    Grid = QTableWidget(nCells.value(), nCells.value())
    OperationsLayout = QVBoxLayout()
    ListLayout = QHBoxLayout()
    IncertingLayout = QVBoxLayout()
    resultLabel = QLabel('Result')
    resultText = QSpinBox()
    resultText.setMaximum(10000)
    operatorLabel = QLabel('Operation')
    operatorText = QLineEdit()
    def AgregarOperaciones():
      listaRes.addItem(str(resultText.value()))
      listaOp.addItem(str(operatorText.text()))
    IncertButton = QPushButton('Agregar')
    IncertButton.clicked.connect(AgregarOperaciones)
    def GetGridItems(NeededValue):
      arr = []
      for i in range(nCells.value()):
        for j in range(nCells.value()):
          #print(Grid.item(i,j).text())
          if(Grid.item(i,j).text() == NeededValue):
            arr += [i,j]
          #print(arr)
      return arr
    def GuardarGenerado():
      filename = QFileDialog.getSaveFileName(CreateWindow,"Guardar KenKen","./kenken.txt","Text Files (*.txt)","Text Files (*.txt)",QFileDialog.Options())
      print(filename)
      file = open(filename[0], 'w')
      for i in range(listaOp.count()):
        arr = GetGridItems(str(i+1))
        file.write(listaRes.item(i).text() +','+listaRes.item(i).text()+ ',' + str(arr))
      file.close()


    GuardarDatos = QPushButton('GuardarKenken')
    GuardarDatos.clicked.connect(GuardarGenerado)
    IncertingLayout.addWidget(resultLabel)
    IncertingLayout.addWidget(resultText)
    IncertingLayout.addWidget(operatorLabel)
    IncertingLayout.addWidget(operatorText)
    IncertingLayout.addWidget(IncertButton)
    IncertingLayout.addWidget(GuardarDatos)
    #List
    ListLayout.addWidget(listaRes)
    ListLayout.addWidget(listaOp)
    OperationsLayout.addLayout(ListLayout)
    OperationsLayout.addLayout(IncertingLayout)
    GridLayout.addWidget(Grid)
    GridLayout.addLayout(OperationsLayout)
    mainLayout.addLayout(GridLayout)
CreateButton = QPushButton('Create Structure')
CreateButton.clicked.connect(CreateGrid)
layout.addWidget(Label_Ken)
layout.addWidget(nCells)
layout.addWidget(CreateButton)

mainLayout.addLayout(layout)

SolveWindow = QWidget()
SolveLayout = QVBoxLayout()
selectorLayout = QHBoxLayout()
ResultLayout = QVBoxLayout()

file = ''
def SearchFile():
    file = QFileDialog().getOpenFileName()
    fileRute.setText(file[0])
def StartSolver():
  MainTab.setGeometry(30,30,1000,600)
  valor = CSP_kenkenSolver(fileRute.text())
  n = len(valor)
  resultTable = QTableWidget(n,n)
  for i in range(n):
    for j in range(n):
      resultTable.setItem(i,j,QTableWidgetItem(str(valor[i][j])))
  print(valor)
  ResultLayout.addWidget(resultTable)
  

## Solve Window
fileSelector = QPushButton('Seleccionar Archivo')
fileSelector.clicked.connect(SearchFile)
fileRute = QLineEdit()
fileRute.setReadOnly(True)
solveButton = QPushButton('Resolver')
solveButton.clicked.connect(StartSolver)
selectorLayout.addWidget(fileSelector)
selectorLayout.addWidget(fileRute)
selectorLayout.addWidget(solveButton) 

SolveLayout.addLayout(selectorLayout)
SolveLayout.addLayout(ResultLayout)

SolveWindow.setLayout(SolveLayout)

CreateWindow.setLayout(mainLayout)

MainTab.addTab(CreateWindow,'Create')
MainTab.addTab(SolveWindow, 'Solve')
MainTab.show()
app.exec_()