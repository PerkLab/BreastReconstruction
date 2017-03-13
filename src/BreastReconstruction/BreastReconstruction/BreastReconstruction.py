import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

#
# BreastReconstruction
#

class BreastReconstruction(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BreastReconstruction" # TODO make this more human readable by adding spaces
    self.parent.categories = ["BreastSurgery"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

#
# BreastReconstructionWidget
#

class BreastReconstructionWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

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
    
    #input model slector

    self.inputModelSelector = slicer.qMRMLNodeComboBox()
    self.inputModelSelector.nodeTypes = [ "vtkMRMLModelNode" ]
    self.inputModelSelector.selectNodeUponCreation = True
    self.inputModelSelector.addEnabled = False
    self.inputModelSelector.removeEnabled = False
    self.inputModelSelector.noneEnabled = False
    self.inputModelSelector.showHidden = False
    self.inputModelSelector.showChildNodeTypes = False
    self.inputModelSelector.setMRMLScene( slicer.mrmlScene )
    self.inputModelSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Model: ", self.inputModelSelector)

    #input Fiducal slector 

    self.inputFiducialSelector = slicer.qSlicerSimpleMarkupsWidget()
    self.inputFiducialSelector.tableWidget().hide()
    self.inputFiducialSelector.setMRMLScene(slicer.mrmlScene)
    self.inputFiducialSelector.setToolTip( "Pick the fiducials to define the region of interest." )
    self.inputFiducialSelector.setNodeBaseName ("BreastFiducials")
    parametersFormLayout.addRow("Input fiducials: ", self.inputFiducialSelector)
    # Enable place multiple marukps by default
    placeWidget = self.inputFiducialSelector.markupsPlaceWidget()
    placeWidget.placeMultipleMarkups = slicer.qSlicerMarkupsPlaceWidget.ForcePlaceMultipleMarkups
    placeWidget.placeModeEnabled = False
    placeWidget.placeModeEnabled = True
    

    #output Model selector
    
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLModelNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Model: ", self.outputSelector)

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
    self.inputFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputModelSelector.currentNode()

  def onApplyButton(self):
    logic = BreastReconstructionLogic()
    logic.run(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode(), self.outputSelector.currentNode())

#
# BreastReconstructionLogic
#

class BreastReconstructionLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def closeModel(self, modelNode, fidList):
    '''
    in this code need to create a rectangle.
    find the points on the edge of the model
    connect edge points to model
    return new model that is now closed
    '''
   
    #create vtkpolydata from input model
    #model = modelNode.GetPolyData()   
    
    
    # VTK connectivity filter
    '''
    surface = modelNode.GetPolyDataConnection()
    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetExtractionModeToLargestRegion()
    connectivityFilter.SetInputConnection(surface)
    surface = connectivityFilter.GetOutputPort()
    modelNode.SetPolyDataConnection(surface)
    modelNode.Update()
    '''

    #Computing normal of first 3 points
    point1 = [0,0,0]
    fidList.GetNthFiducialPosition(0,point1)
    point2 = [0,0,0]
    fidList.GetNthFiducialPosition(1,point2)
    point3 = [0,0,0]
    fidList.GetNthFiducialPosition(2,point3)
    normal = [0,0,0]
    plane = vtk.vtkTriangle()
    plane.ComputeNormal(point1,point2,point3,normal)
    normal[0] = normal[0] * 1
    normal[1] = normal[1] * 1
    normal[2] = normal[2] * 1
    
    #create extrusion fliter to close model
    model = modelNode.GetPolyData() 
    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.SetInputData(model)
    extrude.SetScaleFactor(100)
    extrude.SetExtrusionTypeToVectorExtrusion()
    #extrude.SetExtrusionTypeToNormalExtrusion()
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

    #Create volume for scene 
    imageSize=[512, 512, 512]
    imageSpacing=[1.0, 1.0, 1.0]
    voxelType=vtk.VTK_UNSIGNED_CHAR
    # Create an empty image volume
    imageData=vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 1)
    thresholder=vtk.vtkImageThreshold()
    thresholder.SetInputData(imageData)
    thresholder.SetInValue(0)
    thresholder.SetOutValue(0)
    # Create volume node
    volumeNode=slicer.vtkMRMLScalarVolumeNode()
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetImageDataConnection(thresholder.GetOutputPort())
    # Add volume to scene
    slicer.mrmlScene.AddNode(volumeNode)
    displayNode=slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()
    volumeNode.SetName("Volume Node")


  def run(self, inputModel, fidList, outMode):
    """
    Run the actual algorithm
    """
    self.closeModel(inputModel, fidList)
    
    logging.info('Processing completed')

    return True


class BreastReconstructionTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_BreastReconstruction1()

  def test_BreastReconstruction1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    self.delayDisplay('Test passed!')
