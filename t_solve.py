
import sys,os
from PySide6.QtWidgets import(
     QApplication, 
     QLabel, 
     QPushButton, 
     QWidget, 
     QLineEdit, 
     QVBoxLayout, 
     QTreeWidget,
     QTreeWidgetItem,
     QMainWindow,
     QSpinBox,
     QComboBox,
     QGridLayout,
     QSpacerItem,
     QSizePolicy
     )
sys.path.append("C:/Users/125715/python")

from fluidlib import thermal_solver as t
from PySide6.QtGui import QColor, QFont
from PySide6.QtCore import Slot,  Qt, Signal
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.image import imread
import numpy as np


class MplCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
            fig, self.ax = plt.subplots(figsize=(width, height), dpi=dpi)
            super(MplCanvas, self).__init__(fig) # Set the parent widget to make it behave like a QWidget

            # This ensures that the canvas expands to fit the widget dimensions
            # self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        
class MainWindow(QMainWindow): 

    def __init__(self): 

        super().__init__()
        self.title = 'Thermal Solver'
        self.setWindowTitle(self.title)
        self.setGeometry(100,25,1200,800)

        
        self.table_widget = thermalSolver()
        self.setCentralWidget(self.table_widget)

        self.show()
class thermalSolver(QWidget): 
    
    def __init__(self):

        super().__init__()
        self.nodes = {}
        self.paths = {}
        self.backendNodes = []
        self.backendPaths = []
        
        self.outerLayout = QGridLayout() #base layout
        self.nodeLayout = QVBoxLayout()
        self.attributeLayout = QVBoxLayout()
        self.pathLayout = QVBoxLayout()
        self.pathAttributeLayout = QVBoxLayout()
        self.timeEvalAndUnitLayout = QVBoxLayout()
        self.plotLayout = QVBoxLayout()
        self.unitLayout = QVBoxLayout()

        #LABELS 
        self.numberOfNodeLabel = QLabel("CURRENT NODE")
        self.temperatureLabel = QLabel("NODE TEMPERATURE, K" )
        self.mediumTypeLabel = QLabel("MEDIUM TYPE")
        self.mediumLabel = QLabel("MEDIUM")
        self.pressureLabel= QLabel("NODE PRESSURE, Pa") #Specify if gage or not
        self.volumeLabel = QLabel("NODE VOLUME, m^3") 
        self.heatGeneratedLabel = QLabel("INTERNAL HEAT GENERATION, W")
        self.emissivityLabel= QLabel("EMISSIVITY")
        self.absorbivityLabel = QLabel("ABSORBTIVITY")
        self.isothermalLabel = QLabel("ISOTHERMAL")
        self.pressureUnitsLabel = QLabel("Pressure Units")
        self.temperatureUnitsLabel = QLabel("Temperature Units")
        self.volumeUnitsLabel = QLabel("Volume Units")
        self.nodeTreeTitle = QLabel("NODE TREE")
        self.pathTreeTitle = QLabel("PATH TREE")
        self.connectedNodesLabel = QLabel("CONNECTED NODES")
        self.nodeA = QLabel("NODE A")
        self.nodeB = QLabel("NODE B")
        self.heatTransferCoefficientLabel = QLabel("HEAT TRANSFER COEFFICIENT, W/m^2-K")
        self.heatTransferAreaLabel = QLabel("HEAT TRANSFER AREA, m^2")
        self.numberofPathsLabel =QLabel("CURRENT PATH")     
        self.timeInitialLabel = QLabel("TIME INITIAL,s")
        self.timeFinalLabel = QLabel("TIME FINAL, s")  
        self.dxLabel = QLabel("LENGTH OF PATH, m")
        self.energyUnitsLabel = QLabel("Energy Units")
        self.areaUnitsLabel = QLabel("Area Units")
        self.distanceUnitsLabel = QLabel("Distance Units")
        self.heatTransferModeLabel = QLabel("HEAT TRANSFER MODE")
        
        #UNIT DROPDOWN MENUS
        self.pressureUnits = QComboBox()
        self.pressureUnits.addItems(["Pa","psig", "psia",  "bar,g","bar,a","atm"])
        self.pressureUnits.currentTextChanged.connect(self.updatePressureUnits)

        self.temperatureUnits = QComboBox()
        self.temperatureUnits.addItems(["K","C","F","R"])
        self.temperatureUnits.currentTextChanged.connect(self.updateTemperatureUnits)
        
        self.volumeUnits = QComboBox()
        self.volumeUnits.addItems(["m^3","mm^3","L","ft^3", "in^3", "gal"])
        self.volumeUnits.currentTextChanged.connect(self.updateVolumeUnits)

        self.energyUnits = QComboBox()
        self.energyUnits.addItems(["W", "Btu"])
        self.energyUnits.currentTextChanged.connect(self.updateEnergyUnits)

        self.areaUnits = QComboBox()
        self.areaUnits.addItems(["m^2", "mm^2", "ft^2","in^2"])
        self.areaUnits.currentTextChanged.connect(self.updateAreaUnits)

        self.distanceUnits = QComboBox()
        self.distanceUnits.addItems(["m", "mm", "ft","in"])
        self.distanceUnits.currentTextChanged.connect(self.updateDistanceUnits)

        

        # ERROR PRINTING 
        self.errorText = QLabel()
        self.errorText.setVisible(False)
        #NODE ATTRIBUTE SELECTIONS
        self.currentNodeSelection = QSpinBox()
        self.currentNodeSelection.setRange(1,1)

        self.temperatureInput = QLineEdit()
        self.temperatureInput.setPlaceholderText("300 ")

        self.mediumTypeSelection = QComboBox()
        self.mediumTypeSelection.addItems(["FLUID","SOLID"])
        self.mediumTypeSelection.setCurrentIndex(-1)
        self.mediumTypeSelection.setPlaceholderText("Select")
        self.mediumTypeSelection.currentTextChanged.connect(self.changeMediumSelectionItems)
        self.mediumTypeSelection.currentTextChanged.connect(self.toggleSolidProperties)
        
        self.mediumSelection = QComboBox()
        self.mediumSelection.setPlaceholderText("Select")
        
        self.pressureInput = QLineEdit()
        self.pressureInput.setPlaceholderText("101325")

        self.volumeInput = QLineEdit()
        self.volumeInput.setPlaceholderText(".001")

        self.heatGeneratedInput = QLineEdit()
        self.heatGeneratedInput.setPlaceholderText("100")

        self.emissivityInput = QLineEdit()
        self.emissivityInput.setPlaceholderText("0.3")

        self.absorbivityInput = QLineEdit()
        self.absorbivityInput.setPlaceholderText("0.2")

        self.isothermalInput = QComboBox()
        self.isothermalInput.addItems(["True", "False"])
        self.isothermalInput.setCurrentIndex(1)

        #PATH ATTRIBUTE SELECTIONS
        self.currentPathSelection = QSpinBox()
        self.currentPathSelection.setRange(1, 20)

        self.heatTransferModeInput = QComboBox()
        self.heatTransferModeInput.addItems(["CONVECTION", "CONDUCTION"])
        self.heatTransferModeInput.setCurrentIndex(-1)
        self.heatTransferModeInput.setPlaceholderText("Select")

        self.heatTransferModeInput.currentTextChanged.connect(self.toggleHeatTransferProperties)
        self.heatTransferCoefficientInput = QLineEdit()
        self.heatTransferCoefficientInput.setPlaceholderText("10")

        self.heatTransferAreaInput = QLineEdit()
        self.heatTransferAreaInput.setPlaceholderText("0.1")

        self.dxInput = QLineEdit()
        self.dxInput.setPlaceholderText(".01")

        #TIME INPUT 
        self.timeInitialInput = QLineEdit()
        self.timeInitialInput.setPlaceholderText("0")

        self.timeFinalInput = QLineEdit()
        self.timeFinalInput.setPlaceholderText("3600")

        #NODE TREE
        self.nodeTree = QTreeWidget()
        self.nodeTree.setColumnCount(3)
        self.nodeTree.setColumnWidth(0, 150)
        self.nodeTree.setHeaderLabels(["ATTRIBUTES", "VALUES", "UNITS"])
        self.nodeLayout.addWidget(self.nodeTreeTitle,alignment=Qt.AlignmentFlag.AlignHCenter)
        self.nodeLayout.addWidget(self.nodeTree)

        #PATH TREE
        self.pathTree = QTreeWidget()
        self.pathTree.setColumnCount(3)
        self.pathTree.setColumnWidth(0,150) 
        self.pathTree.setHeaderLabels(["ATRRIBUTES", "VALUES", "UNITS"])
        self.pathLayout.addWidget(self.pathTreeTitle,alignment=Qt.AlignmentFlag.AlignHCenter)
        self.pathLayout.addWidget(self.pathTree)

        #NODE AND PATH SELECTIONS
        self.connectedNodesSelection1 = QComboBox()
        self.connectedNodesSelection2 = QComboBox()
        self.connectedNodesSelection1.setCurrentIndex(-1)
        self.connectedNodesSelection2.setCurrentIndex(-1) 
        self.connectedNodesSelection1.currentIndexChanged.connect(self.changeNodeSelectionItems)
        
        #SET NODE ATTRIBUTE INPUT LAYOUT
        self.attributeLayout.addWidget(self.numberOfNodeLabel)
        self.attributeLayout.addWidget(self.currentNodeSelection)

        self.attributeLayout.addWidget(self.isothermalLabel)
        self.attributeLayout.addWidget(self.isothermalInput)
        
        self.attributeLayout.addWidget(self.temperatureLabel)
        self.attributeLayout.addWidget(self.temperatureInput)

        self.attributeLayout.addWidget(self.volumeLabel)
        self.attributeLayout.addWidget(self.volumeInput)

        self.attributeLayout.addWidget(self.mediumTypeLabel)
        self.attributeLayout.addWidget(self.mediumTypeSelection)
        self.attributeLayout.addWidget(self.mediumLabel)
        self.attributeLayout.addWidget(self.mediumSelection)
        
        self.attributeLayout.addWidget(self.pressureLabel)
        self.attributeLayout.addWidget(self.pressureInput)
        
        self.attributeLayout.addWidget(self.heatGeneratedLabel)
        self.attributeLayout.addWidget(self.heatGeneratedInput)
        
        self.attributeLayout.addWidget(self.emissivityLabel)
        self.attributeLayout.addWidget(self.emissivityInput)
        
        self.attributeLayout.addWidget(self.absorbivityLabel)
        self.attributeLayout.addWidget(self.absorbivityInput)
        
     
        #SET PATH ATTRIBUTE INPUT LAYOUT
        self.pathAttributeLayout.addWidget(self.numberofPathsLabel)
        self.pathAttributeLayout.addWidget(self.currentPathSelection)

        self.pathAttributeLayout.addWidget(self.heatTransferModeLabel)
        self.pathAttributeLayout.addWidget(self.heatTransferModeInput)

        self.pathAttributeLayout.addWidget(self.heatTransferAreaLabel)
        self.pathAttributeLayout.addWidget(self.heatTransferAreaInput)

        self.pathAttributeLayout.addWidget(self.heatTransferCoefficientLabel) 
        self.pathAttributeLayout.addWidget(self.heatTransferCoefficientInput)

        self.pathAttributeLayout.addWidget(self.dxLabel)
        self.pathAttributeLayout.addWidget(self.dxInput)

        self.pathAttributeLayout.addWidget(self.connectedNodesLabel)     
        self.pathAttributeLayout.addWidget(self.nodeA)
        self.pathAttributeLayout.addWidget(self.connectedNodesSelection1)
        self.pathAttributeLayout.addWidget(self.nodeB)
        self.pathAttributeLayout.addWidget(self.connectedNodesSelection2)

        #SET TIME INPUT LAYOUT
        
        self.timeEvalAndUnitLayout.addWidget(self.pressureUnitsLabel)
        self.timeEvalAndUnitLayout.addWidget(self.pressureUnits)
        self.timeEvalAndUnitLayout.addWidget(self.temperatureUnitsLabel)
        self.timeEvalAndUnitLayout.addWidget(self.temperatureUnits)
        self.timeEvalAndUnitLayout.addWidget(self.volumeUnitsLabel)
        self.timeEvalAndUnitLayout.addWidget(self.volumeUnits)
        self.timeEvalAndUnitLayout.addWidget(self.energyUnitsLabel)
        self.timeEvalAndUnitLayout.addWidget(self.energyUnits)
        self.timeEvalAndUnitLayout.addWidget(self.areaUnitsLabel) 
        self.timeEvalAndUnitLayout.addWidget(self.areaUnits)
        self.timeEvalAndUnitLayout.addWidget(self.distanceUnitsLabel)
        self.timeEvalAndUnitLayout.addWidget(self.distanceUnits)
        
        self.timeEvalAndUnitLayout.addWidget(self.timeInitialLabel)
        self.timeEvalAndUnitLayout.addWidget(self.timeInitialInput)
        self.timeEvalAndUnitLayout.addWidget(self.timeFinalLabel)
        self.timeEvalAndUnitLayout.addWidget(self.timeFinalInput)
        self.timeEvalAndUnitLayout.addWidget(self.errorText)

        #PLOT LAYOUT
        self.canvas = MplCanvas(self,width=5, height=2, dpi = 100)
        self.canvas.ax.set_title("Node Temperature vs Time")
        self.canvas.ax.grid(True)
        self.plotLayout.addWidget(self.canvas)
        
        #BUTTONS
        
        self.updateNodeButton = QPushButton("Update Node")

        self.updateNodeButton.pressed.connect(self.updateNode)
        self.updateNodeButton.pressed.connect(self.setNodeSelection)
        self.attributeLayout.addWidget(self.updateNodeButton)

        self.removeNodeButton = QPushButton("Remove Node")
        self.removeNodeButton.pressed.connect(self.removeNode)
        self.attributeLayout.addWidget(self.removeNodeButton)

        self.updatePathButton = QPushButton("Update Path")
        self.updatePathButton.pressed.connect(self.updatePath)

        self.removePathButton = QPushButton("Remove Path")  
        self.removePathButton.pressed.connect(self.removePath)
    
        self.pathAttributeLayout.addWidget(self.updatePathButton)
        self.pathAttributeLayout.addWidget(self.removePathButton)

        self.solveButton = QPushButton("SOLVE")
        self.timeEvalAndUnitLayout.addWidget(self.solveButton)
        self.solveButton.pressed.connect(self.solve)

        self.clearAllButton = QPushButton("CLEAR ALL")
        self.clearAllButton.pressed.connect(self.clearAll)
        self.attributeLayout.addWidget(self.clearAllButton)

        #SET OUTER LAYOUT
        self.outerLayout.addLayout(self.attributeLayout, 0,0)
        self.outerLayout.addLayout(self.nodeLayout,0,1)
        self.outerLayout.addLayout(self.pathLayout, 1, 1)
        self.outerLayout.addLayout(self.pathAttributeLayout,1,0)
        self.outerLayout.addLayout(self.timeEvalAndUnitLayout, 0,2)
        self.outerLayout.addLayout(self.plotLayout,1,2)
        
        self.setLayout(self.outerLayout)
        

    #Set Node Selection based on foregoing Node Selection for path
    def setNodeSelection(self): 
        self.connectedNodesSelection1.clear()
        top_level_count = self.nodeTree.topLevelItemCount()
        self.connectedNodesSelection1.addItem(self.nodeTree.topLevelItem(0).text(0)) if self.nodeTree.topLevelItemCount() != 0 else ''
        for i in range(top_level_count): 
            node_name = self.nodeTree.topLevelItem(i).text(0)
            if self.connectedNodesSelection1.itemText(i) == node_name or self.nodeTree.topLevelItemCount() == 0:
                pass
            else:
                self.connectedNodesSelection1.addItem(node_name)

    '''Add or update node in node tree'''
    def updateNode(self):
        try:
            items = []
            self.nodeTree.clear()
            attributes = [
                "T," + (self.temperatureInput.text()) +':'+self.temperatureUnits.currentText(), 
                "MEDIUM,"+self.mediumSelection.currentText()+':', 
                "PRESSURE,"+ (self.pressureInput.text()) + ':'+self.pressureUnits.currentText() if self.mediumTypeSelection.currentText() == "FLUID" else '' , 
                "HEAT GENERATED," + (self.heatGeneratedInput.text()) + ':'+self.energyUnits.currentText() if self.heatGeneratedInput.text() != '' else '', 
                "VOLUME," + (self.volumeInput.text()) +':'+ self.volumeUnits.currentText(), 
                "EMISSIVITY," + (self.emissivityInput.text()) +':'if self.mediumTypeSelection.currentText() == "SOLID" else '', 
                "ABSORPTIVITY,"+(self.absorbivityInput.text()) +':'if self.mediumTypeSelection.currentText() == "SOLID" else '',
                "ISOTHERMAL," + (self.isothermalInput.currentText()) + ':'
            ]
            '''Stores node in backend nodes and updates nodes using identifier attribute'''
            self.backendNodes = [node for node in self.backendNodes if node.identifier != self.currentNodeSelection.value()]
            self.backendNodes.append(t.Node(
                T = float(self.temperatureInput.text()) if self.temperatureInput.text() != '' else 300, 
                medium= self.mediumSelection.currentText(),
                medium_type=self.mediumTypeSelection.currentText() if self.mediumTypeSelection.currentText() != '' else "SOLID",
                Pressure=float(self.pressureInput.text()) if self.mediumTypeSelection.currentText() == 'FLUID' else 0.0,
                Eg = float(self.heatGeneratedInput.text()) if self.heatGeneratedInput.text() != '' else 0.0, 
                V= float(self.volumeInput.text()) if self.volumeInput.text() != '' else 0.0, 
                #remove emissivity input
                isothermal= eval(self.isothermalInput.currentText()) if self.isothermalInput.currentText() != '' else False,
                identifier=self.currentNodeSelection.value()
                ))

            #Create Node Dictionary
            self.nodes['Node ' + str(self.currentNodeSelection.value())] = attributes
            self.nodes = dict(sorted(self.nodes.items()))
            print(self.nodes)
            print(self.backendNodes)
            for key,values in self.nodes.items():
                if key == '': 
                    pass
                else:
                    nodeAssignment = QTreeWidgetItem([key])
                for value in values:
                    if value == '':
                        pass
                    else:
                        label = value.split(",")[0]
                        val = value.split(':')[0].split(',')[-1]
                        units = value.split(':')[-1]
                        child = QTreeWidgetItem([label,val,units])
                        nodeAssignment.addChild(child)
                items.append(nodeAssignment)
        
            '''if node exists in an already determined path (i.e in path tree), then update path(s) with current node'''
            for node in self.backendNodes:
                if node.identifier == int(self.currentNodeSelection.text()):
                    updatedNode = node
            for i in range(len(self.backendPaths)):
                node1 = self.backendPaths[i].nodeA
                node2 = self.backendPaths[i].nodeB
                if updatedNode != node1 and updatedNode.identifier == node1.identifier:
                    self.backendPaths[i].nodeA = updatedNode
                elif updatedNode != node2 and updatedNode.identifier == node2.identifier:
                    self.backendPaths[i].nodeB = updatedNode

            self.currentNodeSelection.setMaximum(int(self.currentNodeSelection.maximum())+1)
            self.currentNodeSelection.setValue(self.currentNodeSelection.value() + 1)
            self.nodeTree.insertTopLevelItems(0, items)
        except:
            self.errorText.setVisible(True)
            self.errorText.setText("Something went wrong with node update. Check node inputs.")

    def updatePath(self): 
        '''First checks if path already exists AND is not updating pre-existing path'''
        nume1 = int(self.connectedNodesSelection1.currentText().split(' ')[-1])
        denom1 = int(self.connectedNodesSelection2.currentText().split(' ')[-1])
        
        for i in range(self.pathTree.topLevelItemCount()):
            nume2 = int(self.paths['Path ' + str(i+1)][0].split(' ')[-1])
            denom2 = int(self.paths['Path ' + str(i+1)][1].split(' ')[-1])   
            prod = (nume1/denom1)*(nume2/denom2) 
            if self.currentPathSelection.text() != self.pathTree.topLevelItem(i).text(0).split(' ')[-1]:
                if prod == (nume2/denom2)**2 or prod == 1:
                    print("Path already exists")
                    return 1
        self.pathTree.clear()
        attributes =[ 
            "Node A," + self.connectedNodesSelection1.currentText(),
            "Node B," + self.connectedNodesSelection2.currentText(),
            "h," + self.heatTransferCoefficientInput.text(),
            "Area," + self.heatTransferAreaInput.text(),
            "dx," + self.dxInput.text()

        ]
        self.backendPaths = [path for path in self.backendPaths if path.identifier != self.currentPathSelection.value()]
        self.paths = dict(sorted(self.paths.items()))
        '''Look for connected nodes using identifier'''
        for i in range(len(self.backendNodes)):
            if self.backendNodes[i].identifier == int(self.connectedNodesSelection1.currentText().split(" ")[-1]):
                node1 = self.backendNodes[i]
            elif self.backendNodes[i].identifier == int(self.connectedNodesSelection2.currentText().split(" ")[-1]):
                node2 = self.backendNodes[i]
       
        self.backendPaths.append(t.Path(
            nodeA= node1,
            nodeB= node2,
            Area = float(self.heatTransferAreaInput.text()) if self.heatTransferAreaInput.text() != '' else 0.01,
            h = float(self.heatTransferCoefficientInput.text()) if self.heatTransferCoefficientInput != '' else 0.0, 
            dx = float(self.dxInput.text()),
            identifier= self.currentPathSelection.value()
        ))

        self.paths['Path ' + str(self.currentPathSelection.value())] = attributes
        items = []
        for key,values in self.paths.items():
                if key == '': 
                    pass
                else:
                    pathAssignment = QTreeWidgetItem([key])
                for value in values:
                    if value == '':
                        pass
                    else:
                        label = value.split(",")[0]
                        val = value.split(",")[-1]
                        child = QTreeWidgetItem([label,val])
                        pathAssignment.addChild(child)
                items.append(pathAssignment)
        self.currentPathSelection.setValue(self.currentPathSelection.value() + 1)
        self.pathTree.insertTopLevelItems(0, items)
        print(self.paths)
        print(self.backendPaths)

    #Sets Selection Medium based on medium type selection
    def changeMediumSelectionItems(self): 
        self.mediumSelection.clear()
        if self.mediumTypeSelection.currentText() == "SOLID": 
            self.mediumSelection.addItems(["SS316", "Al6061", "Al7075"])
            self.mediumSelection.setCurrentIndex(-1)

        elif self.mediumTypeSelection.currentText() == "FLUID":
            self.mediumSelection.addItems(["HYDROGEN", "ARGON", "NITROGEN"])
            self.mediumSelection.setCurrentIndex(-1)

    #Change connected node selection items in 2nd connected node box
    def changeNodeSelectionItems(self): 
        self.connectedNodesSelection2.clear()
        for i in range(self.nodeTree.topLevelItemCount()):
            if self.connectedNodesSelection1.currentText() == self.nodeTree.topLevelItem(i).text(0) or self.nodeTree.topLevelItemCount() ==0:
                pass
            else:
                self.connectedNodesSelection2.addItem(self.nodeTree.topLevelItem(i).text(0))

    #Adds/Removes Input Boxes depending on selected medium type
    def toggleSolidProperties(self):
        if self.mediumTypeSelection.currentText() == "SOLID": 
            self.pressureLabel.setVisible(False)
            self.pressureInput.setVisible(False)
            self.emissivityLabel.setVisible(True)
            self.emissivityInput.setVisible(True)
            self.absorbivityLabel.setVisible(True)
            self.absorbivityInput.setVisible(True)

        elif self.mediumTypeSelection.currentText() == "FLUID":
            self.pressureLabel.setVisible(True)
            self.pressureInput.setVisible(True)
            self.emissivityLabel.setVisible(False)
            self.emissivityInput.setVisible(False)
            self.absorbivityLabel.setVisible(False)
            self.absorbivityInput.setVisible(False)
   
    #Solves Problem
    def solve(self): 
        try:
            for node in self.backendNodes:
                node.connectedPaths = []
            self.canvas.ax.clear()
            tspan = [float(self.timeInitialInput.text()), float(self.timeFinalInput.text())]
            teval = np.linspace(tspan[0], tspan[1],1000)
            a = 0
            a = t.T_vs_t(tspan, teval, self.backendPaths, self.backendNodes)
            time = a.t
            y =a.y
            legend = []
            for i in range(len(y[:,0])):
                self.canvas.ax.plot(time, y[i,:])
                legend.append(round(a.y[i,0],1))
            self.canvas.ax.set_xlabel("Time ")
            self.canvas.ax.set_ylabel("Temperature")
            self.canvas.ax.set_title("Node Temperature vs Time")
            self.canvas.ax.grid(True)
            self.canvas.ax.legend(legend)
            self.canvas.draw()
            self.errorText.setVisible(False)
        except:
            self.errorText.setVisible(True)
            self.errorText.setText("Solve Inputs are incorrect/incomplete")
   
    #Removes selected node referencing selectedNode box
    def removeNode(self): 
            for i in range(self.nodeTree.topLevelItemCount()):
                if self.currentNodeSelection.text() == self.nodeTree.topLevelItem(i).text(0).split(" ")[-1]:
                    
                    for j in range(len(self.backendNodes)):
                        if self.backendNodes[j].identifier == int(self.currentNodeSelection.text()): 
                            self.backendNodes.remove(self.backendNodes[j])
                            break
                    
                    item = self.nodeTree.takeTopLevelItem(i); del item
                    self.nodes.pop('Node ' + self.currentNodeSelection.text())
                    self.connectedNodesSelection1.removeItem(int(i) - 1)
                    self.currentNodeSelection.setMaximum(int(self.nodeTree.topLevelItem(self.nodeTree.topLevelItemCount()-1).text(0).split(' ')[-1])+1) if self.nodeTree.topLevelItemCount() != 0 else 1
                    self.currentNodeSelection.setValue(self.currentNodeSelection.value() - 1)
                    self.setNodeSelection()
                    self.changeNodeSelectionItems()
                    print(self.nodes)
                    print(self.backendNodes)
                    break
            else: 
                print("Node does not exist")
        
    def removePath(self): 
            a = self.pathTree.topLevelItem(int(self.currentPathSelection.text()) - 1).text(0).split(' ')[-1]
            if self.currentPathSelection.text() == a: 
                item = self.pathTree.takeTopLevelItem(int(self.currentPathSelection.text()) - 1); del item; 
                self.paths.pop('Path ' + self.currentPathSelection.text())
                self.backendPaths.pop(int(self.currentPathSelection.text())-1)
                self.currentPathSelection.setMaximum(self.currentPathSelection.maximum() - 1)
                self.currentPathSelection.setValue(self.currentPathSelection.value() -1)
            else:
                print("Path does not exist")
            print(self.paths)
            print(self.backendPaths)
    
    def clearAll(self):
        for widget in self.findChildren(QLineEdit):
            widget.clear()
        for widget in self.findChildren(QTreeWidget):
            widget.clear()
        self.backendNodes = []
        self.backendPaths =[]
        self.nodes = {}
        self.paths = {}
        self.errorText.setVisible(False)
        self.canvas.ax.clear()
        self.currentNodeSelection.setRange(1,1)
        self.currentPathSelection.setValue(1)
   
    def updateTemperatureUnits(self):
        self.temperatureLabel.setText(
            "NODE TEMPERATURE, " + f'{self.temperatureUnits.currentText()}'
        )
    def updateVolumeUnits(self):
        self.volumeLabel.setText(
            "NODE VOLUME, " + f'{self.volumeUnits.currentText()}'
        )
    def updatePressureUnits(self):
        self.pressureLabel.setText(
            "NODE PRESSURE, " + f'{self.pressureUnits.currentText()}'
        )
    def updateEnergyUnits(self):
        self.heatGeneratedLabel.setText(
            "INTERNAL HEAT GENERATION, " + f'{self.energyUnits.currentText()}'
        )
    def updateAreaUnits(self):
        self.heatTransferAreaLabel.setText(
            "HEAT TRANSFER AREA, " + f'{self.areaUnits.currentText()}'
        )
    def updateDistanceUnits(self):
        self.dxLabel.setText(
            "LENGTH OF PATH, " + f'{self.distanceUnits.currentText()}'
        )
    def toggleHeatTransferProperties(self):
        if self.heatTransferModeInput.currentText() == "CONDUCTION": 
            self.heatTransferCoefficientLabel.setVisible(False)
            self.heatTransferCoefficientInput.setVisible(False)
            self.dxLabel.setVisible(True)
            self.dxInput.setVisible(True)
        else: 
            self.heatTransferCoefficientLabel.setVisible(True)
            self.heatTransferCoefficientInput.setVisible(True)
            self.dxLabel.setVisible(False)
            self.dxInput.setVisible(False)
if __name__ == '__main__': 
    app = QApplication([])
    thermalSolve = MainWindow()
    thermalSolve.show()
    sys.exit(app.exec())

