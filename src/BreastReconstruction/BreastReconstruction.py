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
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    It performs a simple thresholding on the input volume and optionally captures a screenshot.
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

    '''
    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # threshold value
    #
    self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    self.imageThresholdSliderWidget.singleStep = 0.1
    self.imageThresholdSliderWidget.minimum = -100
    self.imageThresholdSliderWidget.maximum = 100
    self.imageThresholdSliderWidget.value = 0.5
    self.imageThresholdSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    parametersFormLayout.addRow("Image threshold", self.imageThresholdSliderWidget)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    '''
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


    self.inputFiducialSelector = slicer.qSlicerSimpleMarkupsWidget()
   # self.inputFiducialSelector.nodeTypes = [ "vtkMRMLMarkupsFiducialNode" ]
    '''
    self.inputFiducialSelector.selectNodeUponCreation = True
    self.inputFiducialSelector.addEnabled = False
    self.inputFiducialSelector.removeEnabled = False
    self.inputFiducialSelector.noneEnabled = False
    self.inputFiducialSelector.showHidden = False
    self.inputFiducialSelector.showChildNodeTypes = False
    '''
    self.inputFiducialSelector.tableWidget().hide()
    self.inputFiducialSelector.setMRMLScene(slicer.mrmlScene)
    self.inputFiducialSelector.setToolTip( "Pick the fiducials to define the region of interest." )
    parametersFormLayout.addRow("Input fiducials: ", self.inputFiducialSelector)
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
    self.applyButton.enabled = self.inputModelSelector.currentNode() and self.inputFiducialSelector.currentNode()

  def onApplyButton(self):
    logic = BreastReconstructionLogic()
   # enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
    #imageThreshold = self.imageThresholdSliderWidget.value
    logic.closeModel(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode())
    #logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), imageThreshold, enableScreenshotsFlag)

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
  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : imageThreshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('BreastReconstructionTest-Start','MyScreenshot',-1)

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
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = BreastReconstructionLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
