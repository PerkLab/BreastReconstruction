import os
import unittest
import vtk, qt, ctk, slicer
import math
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np

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
    self.parent.contributors = ["Rachael House (Perklab, Queen's University)"] # replace with "Firstname Lastname (Organization)"
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

    if self.inputFiducialSelector.currentNode() == None:
        logging.error('Please enter input fiducials')
    else: 
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
  def FiducialsToPolyData(self, fiducials, polyData):
    #create polydata from fiducial list 
    points = vtk.vtkPoints()
    n = fiducials.GetNumberOfFiducials()
    for fiducialIndex in range( 0, n ):
        p = [0, 0, 0]
        fiducials.GetNthFiducialPosition( fiducialIndex, p )
        points.InsertNextPoint( p )
    tempPolyData = vtk.vtkPolyData()
    tempPolyData.SetPoints( points )

    vertex = vtk.vtkVertexGlyphFilter()
    vertex.SetInputData( tempPolyData )
    vertex.Update()

    polyData.ShallowCopy(vertex.GetOutput())

  def LeastSquaresPlane(self, modelNode, fidList, plane):
    '''
    create a least square fitted plane to the points 
    create plane by finding center of mass and direction vector 
    ''' 
    #Store fiducials as PolyData 
    PointsPolyData = vtk.vtkPolyData()
    self.FiducialsToPolyData(fidList, PointsPolyData)
    NumberOfPoints = PointsPolyData.GetNumberOfPoints()
    #Compute the center of mass of all points 
    CenterOfMass = vtk.vtkCenterOfMass()
    CenterOfMass.SetInputData(PointsPolyData)
    CenterOfMass.SetUseScalarsAsWeights(False)
    CenterOfMass.Update()
    center = CenterOfMass.GetCenter()

    tempPlane = vtk.vtkPlane()
    tempPlane.SetOrigin(center)

    plane.SetOrigin(center)
    #Chose 3 of the 8 points to compute the normal from and compute error
    #choose plane with least errror 
    bestDistance = float('inf')
    for i in range(NumberOfPoints):
        for j in range(1, NumberOfPoints):
            for k in range(2, NumberOfPoints):
                if abs(i - j) > 1 and abs(i- k) > 1 and abs(j-k)> 1 and abs(i - j) < (NumberOfPoints - 1) and abs(i- k) < (NumberOfPoints -1):
                    triangle = vtk.vtkTriangle()
                    pi = PointsPolyData.GetPoint(i)
                    pj = PointsPolyData.GetPoint(j)
                    pk = PointsPolyData.GetPoint(k)
                    normalVector = [0.0, 0.0, 0.0]
                    triangle.ComputeNormal(pi, pj, pk, normalVector)
                    tempPlane.SetNormal(normalVector)
                    distance = 0
                    for p in range(NumberOfPoints):
                        point = PointsPolyData.GetPoint(p)
                        distance = distance + tempPlane.DistanceToPlane(point)
                    averageDistance = distance/NumberOfPoints
                    if averageDistance < bestDistance:
                        bestDistance = averageDistance
                        plane.SetNormal(normalVector)
    print(normalVector)

  def ClosedInputSurface(self, modelNode, fidList):

    #Create closed surface for clipping
    extrude = vtk.vtkLinearExtrusionFilter()
    plane = vtk.vtkPlane()
    self.LeastSquaresPlane(modelNode,fidList, plane)
    normVec = plane.GetNormal()
    InputModel = modelNode.GetPolyData() 
    extrude.SetInputData(InputModel)
    extrude.SetScaleFactor(100)
    extrude.SetExtrusionTypeToVectorExtrusion()
    #extrude.SetExtrusionTypeToNormalExtrusion()
    reversedNorm0 = normVec[0]
    reversedNorm1 = normVec[1]
    reversedNorm2 = normVec[2]
    extrude.SetVector(reversedNorm0, reversedNorm1, reversedNorm2)
    extrude.CappingOn()
    extrude.Update()

    # modelsLogic = slicer.modules.models.logic()
    # Model = modelsLogic.AddModel(extrude.GetOutputPort()) 
    # Model.GetDisplayNode().SetVisibility(True)
    # Model.SetName("ClosedBreast")
    # Model.GetDisplayNode().BackfaceCullingOff()

    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.SetInputConnection(extrude.GetOutputPort())
    featureEdges.Update()
    numberOfOpenEdges = featureEdges.GetOutput().GetNumberOfCells()
 
    if(numberOfOpenEdges > 0):
        print("Surface is not closed (closed surface)")
    
    else:
        print("surface is closed (closed surface)")
 



  def createCroppedModel(self, modelNode, fidList):
     #vtkClipPolyData, using defined least-squares plane to crop model
    modelsLogic = slicer.modules.models.logic()

    plane = vtk.vtkPlane()
    self.LeastSquaresPlane(modelNode, fidList, plane)
    InputModel = modelNode.GetPolyData()
    clippedInput = vtk.vtkClipPolyData()
    clippedInput.SetInputData(InputModel)
    clippedInput.SetClipFunction(plane)
    clippedInput.SetValue(0)
    clippedInput.SetInsideOut(True)
    clippedInput.Update()

    cutterPlane = vtk.vtkCutter()
    cutterPlane.SetCutFunction(plane)
    cutterPlane.SetInputData(clippedInput.GetOutput())
    cutterPlane.Update()

    cutterModel = vtk.vtkPolyData()
    cutterModel = cutterPlane.GetOutput()
    surfPlane = vtk.vtkSurfaceReconstructionFilter()
    surfPlane.SetInputData(cutterModel)
    cfPlane = vtk.vtkContourFilter()
    cfPlane.SetInputConnection(surfPlane.GetOutputPort())
    cfPlane.SetValue(0,0.0)
    reversePlane = vtk.vtkReverseSense()
    reversePlane.SetInputConnection(cfPlane.GetOutputPort())
    reversePlane.ReverseCellsOn()
    reversePlane.ReverseNormalsOn()

    reversePolyData = reversePlane.GetOutput()

    # modelsLogic = slicer.modules.models.logic()
    # PlaneModel = modelsLogic.AddModel(reversePlane.GetOutputPort()) 
    # PlaneModel.GetDisplayNode().SetVisibility(False)
    # PlaneModel.SetName("PlaneModel")
    # PlaneModel.GetDisplayNode().BackfaceCullingOff()

    ######################################################
    PointsPolyData = vtk.vtkPolyData()
    self.FiducialsToPolyData(fidList, PointsPolyData)

    loop = vtk.vtkImplicitSelectionLoop()
    loop.SetLoop(PointsPolyData.GetPoints())

    #Now clip with loop
    clipperLoop = vtk.vtkClipPolyData()
    clipperLoop.SetClipFunction(loop)
    clipperLoop.SetInputData(reversePolyData)
    clipperLoop.SetInsideOut(True)
    clipperLoop.Update()

    #clip input with loop
    clipperLoop2 = vtk.vtkClipPolyData()
    clipperLoop2.SetClipFunction(loop)
    clipperLoop2.SetInputData(clippedInput.GetOutput())
    clipperLoop2.SetInsideOut(True)
    clipperLoop2.Update()

    ###########################################################

    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(clipperLoop2.GetOutput())
    appendFilter.AddInputData(clipperLoop.GetOutput())
    appendFilter.Update()

    # modelsLogic = slicer.modules.models.logic()
    # AppendModel = modelsLogic.AddModel(appendFilter.GetOutputPort()) 
    # AppendModel.GetDisplayNode().SetVisibility(True)
    # AppendModel.SetName("AppendModel")
    # AppendModel.GetDisplayNode().BackfaceCullingOff()

    #############################################################
    rotationExtrude = vtk.vtkLinearExtrusionFilter()
    rotationExtrude.SetInputData(clipperLoop2.GetOutput())
    rotationExtrude.SetScaleFactor(100)
    rotationExtrude.CappingOn()
    normVec = plane.GetNormal()
    rotationExtrude.SetVector(normVec[0], normVec[1], normVec[2])
    rotationExtrude.Update()

    #clip extruded breast with plane
    clipperLoop3 = vtk.vtkClipPolyData()
    clipperLoop3.SetClipFunction(plane)
    clipperLoop3.SetInputData(rotationExtrude.GetOutput())
    clipperLoop3.SetInsideOut(True)
    clipperLoop3.Update()

    boundaryEdges = vtk.vtkFeatureEdges()
    boundaryEdges.SetInputData(clipperLoop3.GetOutput())
    boundaryEdges.BoundaryEdgesOn()
    boundaryEdges.FeatureEdgesOff()
    boundaryEdges.NonManifoldEdgesOff()
    boundaryEdges.ManifoldEdgesOff()

    boundaryStrips = vtk.vtkStripper()
    boundaryStrips.SetInputConnection(boundaryEdges.GetOutputPort())
    boundaryStrips.Update()

    boundaryPoly = vtk.vtkPolyData()
    boundaryPoly.SetPoints(boundaryStrips.GetOutput().GetPoints())
    boundaryPoly.SetPolys(boundaryStrips.GetOutput().GetLines())

    appendFilter2 = vtk.vtkAppendPolyData()
    appendFilter2.AddInputData(boundaryPoly)
    appendFilter2.AddInputData(clipperLoop3.GetOutput())
    appendFilter2.Update()

    cleanFinal = vtk.vtkCleanPolyData()
    cleanFinal.SetInputConnection(appendFilter2.GetOutputPort())
    cleanFinal.Update()

    massProperties = vtk.vtkMassProperties()
    massProperties.SetInputConnection(cleanFinal.GetOutputPort())
    volume = massProperties.GetVolume()
    volume = volume/ 1000
    volume = round(volume,2)
    print('Volume in cc')
    print(volume)
    surfaceArea = massProperties.GetSurfaceArea()
    surfaceArea = surfaceArea / 100
    surfaceArea = round(surfaceArea, 2)
    print('Surface Area in cm^2')
    print(surfaceArea)

    ###################################################################


    finalModel = modelsLogic.AddModel(cleanFinal.GetOutputPort()) 
    finalModel.GetDisplayNode().SetVisibility(True)
    finalModel.SetName("CroppedClosedBreast")
    finalModel.GetDisplayNode().BackfaceCullingOff()


    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.SetInputData(cleanFinal.GetOutput())
    featureEdges.Update()
    numberOfOpenEdges = featureEdges.GetOutput().GetNumberOfCells()
 
    if(numberOfOpenEdges > 0):
        print("Surface is not closed (final)")
        print(numberOfOpenEdges)
    
    else:
        print("surface is closed (final)")


  def AddVolumeNode(self):
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
    self.ClosedInputSurface(inputModel, fidList)
    self.createCroppedModel(inputModel, fidList)
    self.AddVolumeNode()
    
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
