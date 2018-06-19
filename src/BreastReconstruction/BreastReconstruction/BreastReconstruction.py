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

    self.planeCheckedBox = qt.QCheckBox("Plane cut")
    self.planeCheckedBox.setCheckState(True)
    # planeCheckedBox.setToolTip("Flip model along its X axis")
    parametersFormLayout.addWidget(self.planeCheckedBox)

    self.curvedCheckedBox = qt.QCheckBox("Curve cut")
    # planeCheckedBox.setToolTip("Flip model along its X axis")
    parametersFormLayout.addWidget(self.curvedCheckedBox)

    self.LeftBreastButton = qt.QPushButton("Left Breast")
    self.LeftBreastButton.enabled = False
    parametersFormLayout.addRow("Left breast computations", self.LeftBreastButton)
   
    self.RightBreastButton = qt.QPushButton("Right Breast")
    self.RightBreastButton.enabled = False
    parametersFormLayout.addRow("Right breast computations", self.RightBreastButton)

    self.VolumeLabelRight = qt.QLabel()
    self.VolumeLabelRight.setText("Volume in cc")
    parametersFormLayout.addRow("Right Breast Volume: ", self.VolumeLabelRight)

    self.VolumeLabelLeft = qt.QLabel()
    self.VolumeLabelLeft.setText("Volume in cc")
    parametersFormLayout.addRow("Left Breast Volume: ", self.VolumeLabelLeft)

    self.SurfaceAreaLabelRight = qt.QLabel()
    self.SurfaceAreaLabelRight.setText("Surface Area in cm^2")
    parametersFormLayout.addRow("Right Breast Surface Area: ", self.SurfaceAreaLabelRight)

    self.SurfaceAreaLabelLeft = qt.QLabel()
    self.SurfaceAreaLabelLeft.setText("Surface Area in cm^2")
    parametersFormLayout.addRow("Left Breast Surface Area: ", self.SurfaceAreaLabelLeft)

#    connections
    self.RightBreastButton.connect('clicked(bool)', self.onRightBreastButton)
    self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.LeftBreastButton.connect('clicked(bool)', self.onLeftBreastButton )

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.RightBreastButton.enabled = self.inputModelSelector.currentNode()
    self.LeftBreastButton.enabled = self.inputModelSelector.currentNode()

  def onRightBreastButton(self):
    if self.planeCheckedBox.isChecked():
        #Use flat cut plane
        logic = BreastReconstructionPlaneLogic()
    elif self.curvedCheckedBox.isChecked():
        #Use curved cut plane
        logic = BreastReconstructionCurveLogic()
    else:
        logging.error('Please select either plane or curved cut')
    
    if self.inputFiducialSelector.currentNode() == None:
        logging.error('Please enter input fiducials')
    elif self.inputModelSelector == None:
        self.error("Please enter input model")
    else: 
        logic.run(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode(), False,self.VolumeLabelRight, self.SurfaceAreaLabelRight)

  def onLeftBreastButton(self, cutPlane):
    if self.planeCheckedBox.isChecked():
        #Use flat cut plane
        logic = BreastReconstructionPlaneLogic()
    elif self.curvedCheckedBox.isChecked():
        #Use curved cut plane
        logic = BreastReconstructionCurveLogic()
    else:
        logging.error('Please select either plane or curved cut')

    if self.inputFiducialSelector.currentNode() == None:
        logging.error('Please enter input fiducials')
    elif self.inputModelSelector == None:
        self.error("Please enter input model")
    else:
        logic.run(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode(), True, self.VolumeLabelLeft, self.SurfaceAreaLabelLeft)

class BreastReconstructionPlaneLogic(ScriptedLoadableModuleLogic):

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
    #Chose 3 of the points to compute the normal from and compute error
    #choose plane with least average distance to the plane (error)
    bestDistance = float('inf')
    for i in range(NumberOfPoints):
        for j in range(1, NumberOfPoints):
            for k in range(2, NumberOfPoints):
                #to ensure no points immediately beside each other are slected
                if abs(i - j) > 1 and abs(i- k) > 1 and abs(j-k)> 1 and abs(i - j) < (NumberOfPoints - 1) and abs(i- k) < (NumberOfPoints -1):
                    triangle = vtk.vtkTriangle()
                    pi = PointsPolyData.GetPoint(i)
                    pj = PointsPolyData.GetPoint(j)
                    pk = PointsPolyData.GetPoint(k)
                    #compute the normal of the 3 selected points 
                    normalVector = [0.0, 0.0, 0.0]
                    triangle.ComputeNormal(pi, pj, pk, normalVector)
                    tempPlane.SetNormal(normalVector)
                    distance = 0
                    for p in range(NumberOfPoints):
                        #compute the distance from each point to the plane
                        point = PointsPolyData.GetPoint(p)
                        distance = distance + abs(tempPlane.DistanceToPlane(point))
                    averageDistance = distance/NumberOfPoints
                    if averageDistance < bestDistance:
                        #select plane if average distance value is lower
                        bestDistance = averageDistance
                        plane.SetNormal(normalVector)

  def ClosedInputSurface(self, modelNode, fidList):
    #Create closed surface for clipping
    #Closed surface is created using the linearExtrusion filter
    #The sirface is extruded in the direction of the breast normal
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

    #Uncomment to addd the closed model to the scene

    modelsLogic = slicer.modules.models.logic()
    Model = modelsLogic.AddModel(extrude.GetOutputPort()) 
    Model.GetDisplayNode().SetVisibility(True)
    Model.SetName("ClosedBreast")
    Model.GetDisplayNode().BackfaceCullingOff()

  def createCroppedModel(self, modelNode, fidList, LeftBreast, volume, surfaceArea):
    # Check which breast volume is being computed for
    if LeftBreast == True:
      name = "ClosedLeftBreast"
    else:
      name = "ClosedRightBreast"

    modelsLogic = slicer.modules.models.logic()
    InputModel = modelNode.GetPolyData()

    # Clip the input model with the plane defined by input points
    plane = vtk.vtkPlane()
    self.LeastSquaresPlane(modelNode, fidList, plane)
    clippedInput = vtk.vtkClipPolyData()
    clippedInput.SetInputData(InputModel)
    clippedInput.SetClipFunction(plane)
    clippedInput.SetValue(0)
    clippedInput.SetInsideOut(LeftBreast)
    clippedInput.Update()

    # create visual respersentation of the plane to add to the scene
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
    cfPlane.SetValue(0, 0.0)
    reversePlane = vtk.vtkReverseSense()
    reversePlane.SetInputConnection(cfPlane.GetOutputPort())
    reversePlane.ReverseCellsOn()
    reversePlane.ReverseNormalsOn()

    reversePolyData = reversePlane.GetOutput()

    modelsLogic = slicer.modules.models.logic()
    PlaneModel = modelsLogic.AddModel(reversePlane.GetOutputPort())
    PlaneModel.GetDisplayNode().SetVisibility(False)
    PlaneModel.SetName("PlaneModel")
    PlaneModel.GetDisplayNode().BackfaceCullingOff()

    # create a loop denfined by the input points
    PointsPolyData = vtk.vtkPolyData()
    self.FiducialsToPolyData(fidList, PointsPolyData)

    loop = vtk.vtkImplicitSelectionLoop()
    loop.SetLoop(PointsPolyData.GetPoints())
    loop.SetNormal(plane.GetNormal())


    # Clip the clipped input model with the loop
    clippedInputWithLoop = vtk.vtkClipPolyData()
    clippedInputWithLoop.SetClipFunction(loop)  # should be loop
    clippedInputWithLoop.SetInputData(clippedInput.GetOutput())
    clippedInputWithLoop.SetInsideOut(True)
    clippedInputWithLoop.Update()


    # close the clippedInputWith loop by using the linearExtrusion filter
    extrudeInputWithLoop = vtk.vtkLinearExtrusionFilter()
    extrudeInputWithLoop.SetInputData(clippedInputWithLoop.GetOutput())
    extrudeInputWithLoop.SetScaleFactor(100)
    extrudeInputWithLoop.CappingOn()
    normVec = plane.GetNormal()
    if LeftBreast == True:
      extrudeInputWithLoop.SetVector((normVec[0]), (normVec[1]), (normVec[2]))
    else:
      extrudeInputWithLoop.SetVector((normVec[0] * -1), (normVec[1] * -1), (normVec[2] * -1))
    extrudeInputWithLoop.Update()


    # Cpmpute Point Normals
    extrudeNormals = vtk.vtkPolyDataNormals()
    extrudeNormals.SetInputData(extrudeInputWithLoop.GetOutput())
    extrudeNormals.ComputePointNormalsOn()
    #Do not use autoorient normals here, it will cause some surfaces not to be closed when
    #clipped with plane
    extrudeNormals.Update()

    if LeftBreast == True:
      plane.SetNormal((normVec[0] * -1), (normVec[1] * -1), (normVec[2] * -1))

    planeCollection = vtk.vtkPlaneCollection()
    planeCollection.AddItem(plane)

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(extrudeNormals.GetOutput())
    clean.Update()

    clipClosedBreast = vtk.vtkClipClosedSurface()
    clipClosedBreast.SetInputConnection(clean.GetOutputPort())
    clipClosedBreast.SetClippingPlanes(planeCollection)
    clipClosedBreast.TriangulationErrorDisplayOn()
    clipClosedBreast.Update()


    # extract the volume and surface area poperties from the closed breast
    massProperties = vtk.vtkMassProperties()
    massProperties.SetInputConnection(clipClosedBreast.GetOutputPort())
    volumeMP = massProperties.GetVolume()
    volumeMP = volumeMP / 1000
    volumeMP = round(volumeMP, 2)
    volume.setText(volumeMP)
    # print('Volume in cc')
    # print(volumeMP)
    surfaceAreaMP = massProperties.GetSurfaceArea()
    surfaceAreaMP = surfaceAreaMP / 100
    surfaceAreaMP = round(surfaceAreaMP, 2)
    surfaceArea.setText(surfaceAreaMP)
    # print('Surface Area in cm^2')
    # print(surfaceAreaMP)

    # add closed breast to the scene
    finalModel = modelsLogic.AddModel(clipClosedBreast.GetOutputPort())
    finalModel.GetDisplayNode().SetVisibility(True)
    finalModel.SetName(name)
    finalModel.GetDisplayNode().BackfaceCullingOff()

    #Ensure the final model is closed so the volume computation is correct
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.SetInputData(clipClosedBreast.GetOutput())
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

  def run(self, inputModel, fidList, LeftBreast, volume, surfaceArea):
    """
    Run the actual algorithm
    """
    self.createCroppedModel(inputModel, fidList, LeftBreast, volume, surfaceArea)
    self.AddVolumeNode()
    
    logging.info('Processing completed')

    return True

class BreastReconstructionCurveLogic(ScriptedLoadableModuleLogic):

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
    #Chose 3 of the points to compute the normal from and compute error
    #choose plane with least average distance to the plane (error)
    bestDistance = float('inf')
    for i in range(NumberOfPoints):
        for j in range(1, NumberOfPoints):
            for k in range(2, NumberOfPoints):
                #to ensure no points immediately beside each other are slected
                if abs(i - j) > 1 and abs(i- k) > 1 and abs(j-k)> 1 and abs(i - j) < (NumberOfPoints - 1) and abs(i- k) < (NumberOfPoints -1):
                    triangle = vtk.vtkTriangle()
                    pi = PointsPolyData.GetPoint(i)
                    pj = PointsPolyData.GetPoint(j)
                    pk = PointsPolyData.GetPoint(k)
                    #compute the normal of the 3 selected points 
                    normalVector = [0.0, 0.0, 0.0]
                    triangle.ComputeNormal(pi, pj, pk, normalVector)
                    tempPlane.SetNormal(normalVector)
                    distance = 0
                    for p in range(NumberOfPoints):
                        #compute the distance from each point to the plane
                        point = PointsPolyData.GetPoint(p)
                        distance = distance + abs(tempPlane.DistanceToPlane(point))
                    averageDistance = distance/NumberOfPoints
                    if averageDistance < bestDistance:
                        #select plane if average distance value is lower
                        bestDistance = averageDistance
                        plane.SetNormal(normalVector)

  def ClosedInputSurface(self, modelNode, fidList):
    #Create closed surface for clipping
    #Closed surface is created using the linearExtrusion filter
    #The sirface is extruded in the direction of the breast normal
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

    #Uncomment to addd the closed model to the scene

    modelsLogic = slicer.modules.models.logic()
    Model = modelsLogic.AddModel(extrude.GetOutputPort()) 
    Model.GetDisplayNode().SetVisibility(True)
    Model.SetName("ClosedBreast")
    Model.GetDisplayNode().BackfaceCullingOff()

  def createCroppedModel(self, modelNode, fidList, LeftBreast, volume, surfaceArea):
     #Check which breast volume is being computed for 
    if LeftBreast == True:
        name = "ClosedLeftBreast"
    else:
        name = "ClosedRightBreast"

    modelsLogic = slicer.modules.models.logic()

    #Clip the input model with the plane defined by input points
    plane = vtk.vtkPlane()
    self.LeastSquaresPlane(modelNode, fidList, plane)
    InputModel = modelNode.GetPolyData()
    clippedInput = vtk.vtkClipPolyData()
    clippedInput.SetInputData(InputModel)
    clippedInput.SetClipFunction(plane)
    clippedInput.SetValue(0)
    clippedInput.SetInsideOut(LeftBreast)
    clippedInput.Update()

    #create visual respersentation of the plane to add to the scene
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

    modelsLogic = slicer.modules.models.logic()
    PlaneModel = modelsLogic.AddModel(reversePlane.GetOutputPort()) 
    PlaneModel.GetDisplayNode().SetVisibility(False)
    PlaneModel.SetName("PlaneModel")
    PlaneModel.GetDisplayNode().BackfaceCullingOff()

    #create a loop denfined by the input points
    PointsPolyData = vtk.vtkPolyData()
    self.FiducialsToPolyData(fidList, PointsPolyData)


    # #****************************************************************************
    #spline model to model back of the chest wall 
    spline = vtk.vtkParametricSpline()
    splinePoints = vtk.vtkPoints()
    for i in range(PointsPolyData.GetNumberOfPoints()):
        splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))
    #should the orgin of the plane by included in the point set or not? 
    #splinePoints.InsertNextPoint(plane.GetOrigin())

    spline.SetPoints(splinePoints)
    spline.ClosedOn()
    functionSource = vtk.vtkParametricFunctionSource()
    functionSource.SetParametricFunction(spline)
    implictSpline =  vtk.vtkImplicitPolyDataDistance()
    implictSpline.SetInput(functionSource.GetOutput())

    #Create a new dataset which is the contour points projected onto the plane
    projectedPoints = vtk.vtkPoints()
    projectedPointsPolyData = vtk.vtkPolyData()
    NumberOfPoints = splinePoints.GetNumberOfPoints()
    for i in range(NumberOfPoints):
        p = splinePoints.GetPoint(i)
        pProj = [0,0,0]
        plane.ProjectPoint(p,pProj)
        projectedPoints.InsertNextPoint(pProj)
    projectedPointsPolyData.SetPoints(projectedPoints)

    spline2 = vtk.vtkParametricSpline()
    spline2.SetPoints(projectedPointsPolyData.GetPoints())
    spline2.ClosedOn()
    functionSource2 = vtk.vtkParametricFunctionSource()
    functionSource2.SetParametricFunction(spline2)

    splineTransform = vtk.vtkThinPlateSplineTransform()
    splineTransform.SetSourceLandmarks(projectedPoints)
    splineTransform.SetTargetLandmarks(PointsPolyData.GetPoints())
    splineTransform.SetBasisToR()

    TransformedPlane = vtk.vtkTransformPolyDataFilter()
    TransformedPlane.SetInputData(reversePolyData)
    TransformedPlane.SetTransform(splineTransform)

    finalModel = modelsLogic.AddModel(TransformedPlane.GetOutputPort()) 
    finalModel.GetDisplayNode().SetVisibility(False)
    finalModel.SetName("transformedPlane")
    finalModel.GetDisplayNode().BackfaceCullingOff()

    implictSplinePlane =  vtk.vtkImplicitPolyDataDistance()
    implictSplinePlane.SetInput(TransformedPlane.GetOutput())

    clippedTransformedPlane = vtk.vtkClipPolyData()
    clippedTransformedPlane.SetClipFunction(implictSplinePlane) #should be loop
    clippedTransformedPlane.SetInputData(InputModel)
    clippedTransformedPlane.SetInsideOut(False)
    clippedTransformedPlane.Update()

    ###NOTE: clipping with the selection loop does not work as it will create an open surface... try to clip with the spline
    loop = vtk.vtkImplicitSelectionLoop()
    loop.SetLoop(PointsPolyData.GetPoints())
    loop.SetNormal(plane.GetNormal())

    #Clip the clipped input model with the loop
    clippedInputWithLoop = vtk.vtkClipPolyData()
    clippedInputWithLoop.SetClipFunction(loop) #should be loop
    clippedInputWithLoop.SetInputData(clippedTransformedPlane.GetOutput())
    clippedInputWithLoop.SetInsideOut(True)
    clippedInputWithLoop.Update()


    #These normals must be flipped so that the back of the breast is orriented the correct way
    Normals = vtk.vtkPolyDataNormals()
    Normals.SetInputData(clippedInputWithLoop.GetOutput())
    Normals.ComputePointNormalsOn()
    Normals.ConsistencyOn()
    Normals.FlipNormalsOn()
    Normals.Update()

    boundaryEdges = vtk.vtkFeatureEdges()
    boundaryEdges.SetInputConnection(Normals.GetOutputPort())
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


    #####*****Find out why this turns back black an why surface is not closed 
    Normals2 = vtk.vtkPolyDataNormals()
    Normals2.SetInputData(boundaryPoly)
    Normals2.ConsistencyOn()
    Normals2.SplittingOn()   

    Normals2.Update()

    appendFilter2 = vtk.vtkAppendPolyData()
    appendFilter2.AddInputData(clippedInputWithLoop.GetOutput())
    appendFilter2.AddInputData(Normals2.GetOutput())
    appendFilter2.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(appendFilter2.GetOutput())
    clean.Update()

    finalModel = modelsLogic.AddModel(clean.GetOutput()) 
    finalModel.GetDisplayNode().SetVisibility(True)
    finalModel.SetName("ClosedModelTransformed")
    finalModel.GetDisplayNode().BackfaceCullingOff()

    print("New Volume and SurfaceA")
    massProperties = vtk.vtkMassProperties()
    massProperties.SetInputData(clean.GetOutput())
    volumeMP = massProperties.GetVolume()
    volumeMP = volumeMP/ 1000
    volumeMP = round(volumeMP,2)
    volume.setText(volumeMP)
    print('Volume in cc')
    print(volumeMP)
    surfaceAreaMP = massProperties.GetSurfaceArea()
    surfaceAreaMP = surfaceAreaMP / 100
    surfaceAreaMP = round(surfaceAreaMP, 2)
    surfaceArea.setText(surfaceAreaMP)
    print('Surface Area in cm^2')
    print(surfaceAreaMP)

    #Ensure the final model is closed so the volume computation is correct
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.SetInputData(clean.GetOutput())
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

  def run(self, inputModel, fidList, LeftBreast, volume, surfaceArea):
    """
    Run the actual algorithm
    """
    self.createCroppedModel(inputModel, fidList, LeftBreast, volume, surfaceArea)
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