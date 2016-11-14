import random
import os
import unittest
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

#
# CloseModel
#

class CloseModel(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "CloseModel" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = [] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# CloseModelWidget
#

class CloseModelWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
   # Get module name by stripping 'Widget' from the class name
  def __init__(self, parent = None):
    self.moduleName = self.__class__.__name__
    if self.moduleName.endswith('Widget'):
      self.moduleName = self.moduleName[:-6]
    settings = qt.QSettings()
    self.developerMode = settings.value('Developer/DeveloperMode').lower() == 'true'
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

    def setupDeveloperSection(self):
      if not self.developerMode:
        return

 

   

  def setup(self):
    # Instantiate and connect default widgets ...
    self.setupDeveloperSection()

  def cleanup(self):
    pass

  def onReload(self):
    """
    ModuleWizard will substitute correct default moduleName.
    Generic reload method for any scripted module.
    """
    slicer.util.reloadScriptedModule(self.moduleName)

  def onReloadAndTest(self):
    try:
      self.onReload()
      test = slicer.selfTests[self.moduleName]
      test()
    except Exception, e:
      import traceback
      traceback.print_exc()
      errorMessage = "Reload and Test: Exception!\n\n" + str(e) + "\n\nSee Python Console for Stack Trace"
      slicer.util.errorDisplay(errorMessage)

  def onEditSource(self):
    filePath = slicer.util.modulePath(self.moduleName)
    qt.QDesktopServices.openUrl(qt.QUrl("file:///"+filePath, qt.QUrl.TolerantMode))

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input model selector
    #
    self.inputModelSelector = slicer.qMRMLNodeComboBox()
    self.inputModelSelector.nodeTypes = ( ("vtkMRMLModelNode"), "" )
    self.inputModelSelector.selectNodeUponCreation = True
    self.inputModelSelector.addEnabled = False
    self.inputModelSelector.removeEnabled = False
    self.inputModelSelector.noneEnabled = False
    self.inputModelSelector.showHidden = False
    self.inputModelSelector.showChildNodeTypes = False
    self.inputModelSelector.setMRMLScene( slicer.mrmlScene )
    self.inputModelSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Model: ", self.inputModelSelector)


    self.inputFiducialSelector = slicer.qMRMLNodeComboBox()
    self.inputFiducialSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.inputFiducialSelector.selectNodeUponCreation = True
    self.inputFiducialSelector.addEnabled = False
    self.inputFiducialSelector.removeEnabled = False
    self.inputFiducialSelector.noneEnabled = False
    self.inputFiducialSelector.showHidden = False
    self.inputFiducialSelector.showChildNodeTypes = False
    self.inputFiducialSelector.setMRMLScene( slicer.mrmlScene )
    self.inputFiducialSelector.setToolTip( "Pick the fiducials to define the region of interest." )
    parametersFormLayout.addRow("Input fiducials: ", self.inputFiducialSelector)

   


    #
    # output model selector
    #
    ''''
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLModelNode"), "" )
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Model: ", self.outputSelector)
    '''

    # Parameters Area
    #
   #
    # Reload and Test area
    # Used during development, but hidden when delivering
    # developer mode is turned off.
    self.reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    self.reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(self.reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(self.reloadCollapsibleButton)

    # reload button
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "ScriptedLoadableModuleTemplate Reload"
    self.reloadButton.connect('clicked()', self.onReload)

    # reload and test button
    self.reloadAndTestButton = qt.QPushButton("Reload and Test")
    self.reloadAndTestButton.toolTip = "Reload this module and then run the self tests."
    self.reloadAndTestButton.connect('clicked()', self.onReloadAndTest)

    # edit python source code
    self.editSourceButton = qt.QPushButton("Edit")
    self.editSourceButton.toolTip = "Edit the module's source code."
    self.editSourceButton.connect('clicked()', self.onEditSource)

    # restart Slicer button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.restartButton = qt.QPushButton("Restart Slicer")
    self.restartButton.toolTip = "Restart Slicer"
    self.restartButton.name = "ScriptedLoadableModuleTemplate Restart"
    self.restartButton.connect('clicked()', slicer.app.restart)

    def createHLayout(elements):
      widget = qt.QWidget()
      rowLayout = qt.QHBoxLayout()
      widget.setLayout(rowLayout)
      for element in elements:
        rowLayout.addWidget(element)
      return widget

    #reloadFormLayout.addWidget(createHLayout([self.reloadButton, self.reloadAndTestButton, self.editSourceButton, self.restartButton]))
    #removed reload and test button for now as there is no test
    reloadFormLayout.addWidget(createHLayout([self.reloadButton, self.editSourceButton, self.restartButton]))

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
   # self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputModelSelector.currentNode() # and self.outputSelector.currentNode()

  def onApplyButton(self):
    logic = CloseModelLogic()
    logic.VolumeCalculation(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode())
#
# CloseModelLogic
#

class CloseModelLogic(ScriptedLoadableModuleLogic):
 def VolumeCalculation(self, modelNode, fidList):
  '''
  in this code need to create a rectangle.
  find the points on the edge of the model
  connect edge points to model
  return new model that is now closed
  
  pointIdList = vtk.vtkIdList()
  edge = vtk.vtkCell()
 
  for i in range(0, modelNode.etNumberOfEdges()):
    {
      edge.InsertNextCell(modelNode.GetEdge(i))
      pointIdList.InsertNextId(edge.GetPointIds())
    }

  cleaner = vtk.vtkCleanPolyData()
  cleaner.SetInput(modelNode)
  cleaner.Update()

  triangleFilter = vtk.vtkTriangleFilter()
  triangleFilter.SetInput(cleaner.GetOutput())
  triangleFilter.Update()

  massProps = vtk.vtkMassProperties()
  massProps.SetInput(triangleFilter.GetOutput())
  massProps.Update()

  modelNode.SurfaceArea = massProps.GetSurfaceArea()
  modelNode.Volume = massProps.GetVolume()
  modelNode.ShapeIndex = massProps.GetNormalizedShapeIndex()  
  '''

  #create vtkpolydata from input model
  model = modelNode.GetPolyData()   

  normals = vtk.vtkPolyDataNormals()
  normals.SetInputData(model)
  normals.ComputePointNormalsOn()
  normals.ComputeCellNormalsOff()
  normals.Update()

  #Computing normal of first 3 points
  plane = vtk.vtkTriangle()
  point1 = [0,0,0]
  fidList.GetNthFiducialPosition(1,point1)
  point2 = [0,0,0]
  fidList.GetNthFiducialPosition(2,point2)
  point3 = [0,0,0]
  fidList.GetNthFiducialPosition(3,point3)
  cubeCenter = [0,0,0]
  fidList.GetNthFiducialPosition(4,cubeCenter)


  normal = [0,0,0]
  plane.ComputeNormal(point1,point2,point3,normal)
  normal[0] = normal[0] * -1
  normal[1] = normal[1] * -1
  normal[2] = normal[2] * -1
  
  #create extrusion fliter to close model
  extrude = vtk.vtkLinearExtrusionFilter()
  extrude.SetInputData(model)
  extrude.SetScaleFactor(250)
  extrude.SetExtrusionTypeToNormalExtrusion()
  extrude.SetVector(normal)
  extrude.CappingOn()

  #create new model using excrusion filter
  newmodel = vtk.vtkPolyData()
  newmodel = extrude.GetOutputPort()

 
  #adding the new model to the scene
  modelsLogic = slicer.modules.models.logic()
  outPutModel = modelsLogic.AddModel(newmodel)
  outPutModel.SetName("Closed Model")
  outPutModel.GetDisplayNode().SetVisibility(True)
  outPutModel.GetDisplayNode().SetColor(0.9, 0.3, 0.9)
  outPutModel.GetDisplayNode().SetOpacity(1.0)  # can play with these settings
  outPutModel.GetDisplayNode().SetAmbient(0.1)
  outPutModel.GetDisplayNode().SetDiffuse(0.9)
  outPutModel.GetDisplayNode().SetSpecular(0.1)
  outPutModel.GetDisplayNode().SetPower(10)
  outPutModel.GetDisplayNode().BackfaceCullingOff()

'''
  box = vtk.vtkCubeSource()
  box.SetCenter(cubeCenter)
  box.SetXLength(100)
  box.SetYLength(100)
  box.SetZLength(100)
  polyTransformToProbe = vtk.vtkTransformPolyDataFilter()
  polyTransformToProbe.SetInputConnection(box.GetOutputPort())

  BoxModel = modelsLogic.AddModel(polyTransformToProbe.GetOutputPort())
  BoxModel.SetName("Box Model")
  BoxModel.GetDisplayNode().SetVisibility(True)
  BoxModel.GetDisplayNode().SetColor(0.9, 0.7, 0.4)
  BoxModel.GetDisplayNode().SetOpacity(0.5)
  BoxModel = slicer.util.getNode("Box Model")

'''

class CloseModelTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  

