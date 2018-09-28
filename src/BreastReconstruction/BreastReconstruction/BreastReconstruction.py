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

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "BreastReconstruction" # TODO make this more human readable by adding spaces
    self.parent.categories = ["BreastSurgery"]
    self.parent.dependencies = []
    self.parent.contributors = ["Rachael House (Perklab, Queen's University)"]
    self.parent.helpText = """
    This module takes a 3D surface scans of a woman's chest and points which outline the breast. 
    The breast is segmented for the rest of the surface scan and the breast volume is computed. 
    """
    self.parent.acknowledgementText = """
    This module was developed in the Perklab at Queen's University.
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
    
    #input model selector
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

    #input Fiducal selector
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
    parametersFormLayout.addWidget(self.planeCheckedBox)

    self.curvedCheckedBox = qt.QCheckBox("Curve cut")
    self.curvedCheckedBox.setCheckState(True)
    parametersFormLayout.addWidget(self.curvedCheckedBox)

    #output model selector
    #output currently not working
    self.outputModelSelector = slicer.qMRMLNodeComboBox()
    self.outputModelSelector.nodeTypes = ["vtkMRMLModelNode"]
    self.outputModelSelector.selectNodeUponCreation = True
    self.outputModelSelector.addEnabled = True
    self.outputModelSelector.renameEnabled = True
    self.outputModelSelector.removeEnabled = True
    self.outputModelSelector.noneEnabled = False
    self.outputModelSelector.showHidden = False
    self.outputModelSelector.showChildNodeTypes = False
    self.outputModelSelector.baseName = "outputBreastModel"
    self.outputModelSelector.setMRMLScene( slicer.mrmlScene )
    parametersFormLayout.addRow("Output Model: ", self.outputModelSelector)

    #Add buttons
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
    self.LeftBreastButton.connect('clicked(bool)', self.onLeftBreastButton)
    self.outputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

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
        planeType = True
    elif self.curvedCheckedBox.isChecked():
        #Use curved cut plane
        planeType = False 
    else:
        logging.error('Please select either plane or curved cut')
    
    if self.inputFiducialSelector.currentNode() == None:
        logging.error('Please enter input fiducials')
    elif self.inputModelSelector == None:
        self.error("Please enter input model")
    else:
        logic = BreastReconstructionLogic() 
        logic.run(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode(), False,self.VolumeLabelRight, self.SurfaceAreaLabelRight, self.outputModelSelector.currentNode(), planeType)

  def onLeftBreastButton(self):
    if self.planeCheckedBox.isChecked():
        #Use flat cut plane
        planeType = True
    elif self.curvedCheckedBox.isChecked():
        #Use curved cut plane
        planeType = False 
    else:
        logging.error('Please select either plane or curved cut')

    if self.inputFiducialSelector.currentNode() == None:
        logging.error('Please enter input fiducials')
    elif self.inputModelSelector == None:
        self.error("Please enter input model")
    else:
        logic = BreastReconstructionLogic()
        logic.run(self.inputModelSelector.currentNode(),self.inputFiducialSelector.currentNode(), True, self.VolumeLabelLeft, self.SurfaceAreaLabelLeft, self.outputModelSelector.currentNode(), planeType)

class BreastReconstructionLogic(ScriptedLoadableModuleLogic):

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
    #brute force approach
    bestDistance = float('inf')
    for i in range(NumberOfPoints):
        for j in range(1, NumberOfPoints):
            for k in range(2, NumberOfPoints):
                #to ensure no points immediately beside each other are selected
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


  def cropWithPlane(self, modelNode, fidList, LeftBreast, volume, surfaceArea, output):
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

    # create visual representation of the plane to add to the scene
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


    # create a loop defined by the input points
    PointsPolyData = vtk.vtkPolyData()
    self.FiducialsToPolyData(fidList, PointsPolyData)





    # #****************************************************************************
    # spline model to model back of the chest wall
    spline = vtk.vtkParametricSpline()
    splinePoints = vtk.vtkPoints()
    for i in range(PointsPolyData.GetNumberOfPoints()):
        splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))

    spline.SetPoints(splinePoints)
    spline.ClosedOn()
    functionSource = vtk.vtkParametricFunctionSource()
    functionSource.SetParametricFunction(spline)
    implictSpline = vtk.vtkImplicitPolyDataDistance()
    implictSpline.SetInput(functionSource.GetOutput())

    # Create a new dataset which is the contour points projected onto the plane
    projectedPoints = vtk.vtkPoints()
    projectedPointsPolyData = vtk.vtkPolyData()
    NumberOfPoints = splinePoints.GetNumberOfPoints()
    for i in range(NumberOfPoints):
        p = splinePoints.GetPoint(i)
        pProj = [0, 0, 0]
        plane.ProjectPoint(p, pProj)
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
    TransformedPlane.SetInputConnection(reversePlane.GetOutputPort())
    TransformedPlane.SetTransform(splineTransform)

    finalModel = modelsLogic.AddModel(TransformedPlane.GetOutputPort())
    finalModel.GetDisplayNode().SetVisibility(False)
    finalModel.SetName("transformedPlane")

    implictSplinePlane = vtk.vtkImplicitPolyDataDistance()
    implictSplinePlane.SetInput(TransformedPlane.GetOutput())
#######################################################################################
    loop = vtk.vtkImplicitSelectionLoop()
    loop.SetLoop(PointsPolyData.GetPoints())
    loop.SetNormal(plane.GetNormal())


    # Clip the clipped input model with the loop
    clippedInputWithLoop = vtk.vtkClipPolyData()
    clippedInputWithLoop.SetClipFunction(implictSplinePlane)  # should be loop
    clippedInputWithLoop.SetInputData(clippedInput.GetOutput())
    clippedInputWithLoop.SetInsideOut(False)
    clippedInputWithLoop.Update()

    # No use the vtkPolyDataConnectivityFilter to extract the largest region
    connectedInput = vtk.vtkPolyDataConnectivityFilter()
    connectedInput.SetInputConnection(clippedInputWithLoop.GetOutputPort())
    connectedInput.SetExtractionModeToLargestRegion()
    connectedInput.Update()


    # close the clippedInputWith loop by using the linearExtrusion filter
    extrudeInputWithLoop = vtk.vtkLinearExtrusionFilter()
    extrudeInputWithLoop.SetInputData(connectedInput.GetOutput())
    extrudeInputWithLoop.SetScaleFactor(100)
    extrudeInputWithLoop.CappingOn()
    normVec = plane.GetNormal()
    if LeftBreast == True:
      extrudeInputWithLoop.SetVector((normVec[0]), (normVec[1]), (normVec[2]))
    else:
      extrudeInputWithLoop.SetVector((normVec[0] * -1), (normVec[1] * -1), (normVec[2] * -1))
    extrudeInputWithLoop.Update()


    # Compute Point Normals
    extrudeNormals = vtk.vtkPolyDataNormals()
    extrudeNormals.SetInputData(extrudeInputWithLoop.GetOutput())
    extrudeNormals.ComputePointNormalsOn()
    #Do not use auto orient normals here, it will cause some surfaces not to be closed when
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

    # extract the volume and surface area properties from the closed breast
    massProperties = vtk.vtkMassProperties()
    massProperties.SetInputConnection(clipClosedBreast.GetOutputPort())
    volumeMP = massProperties.GetVolume()
    volumeMP = volumeMP / 1000
    volumeMP = round(volumeMP, 2)
    volume.setText(volumeMP)

    surfaceAreaMP = massProperties.GetSurfaceArea()
    surfaceAreaMP = surfaceAreaMP / 100
    surfaceAreaMP = round(surfaceAreaMP, 2)
    surfaceArea.setText(surfaceAreaMP)


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

  def cropWithCurve(self, modelNode, fidList, LeftBreast, volume, surfaceArea, output):
     #Check which breast volume is being computed for 
     if LeftBreast == True:
       name = "ClosedLeftBreast"
     else:
       name = "ClosedRightBreast"

     modelsLogic = slicer.modules.models.logic()

     # Clip the input model with the plane defined by input points
     plane = vtk.vtkPlane()
     self.LeastSquaresPlane(modelNode, fidList, plane)
     InputModel = modelNode.GetPolyData()
     clippedInput = vtk.vtkClipPolyData()
     clippedInput.SetInputData(InputModel)
     clippedInput.SetClipFunction(plane)
     clippedInput.SetValue(0)
     clippedInput.SetInsideOut(LeftBreast)
     clippedInput.Update()

     # create visual representation of the plane to add to the scene
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

     modelsLogic = slicer.modules.models.logic()

     # create a loop defined by the input points
     PointsPolyData = vtk.vtkPolyData()
     self.FiducialsToPolyData(fidList, PointsPolyData)

     # #****************************************************************************
     # spline model to model back of the chest wall
     spline = vtk.vtkParametricSpline()
     splinePoints = vtk.vtkPoints()
     for i in range(PointsPolyData.GetNumberOfPoints()):
       splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))

     spline.SetPoints(splinePoints)
     spline.ClosedOn()
     functionSource = vtk.vtkParametricFunctionSource()
     functionSource.SetParametricFunction(spline)
     implictSpline = vtk.vtkImplicitPolyDataDistance()
     implictSpline.SetInput(functionSource.GetOutput())

     # Create a new dataset which is the contour points projected onto the plane
     projectedPoints = vtk.vtkPoints()
     projectedPointsPolyData = vtk.vtkPolyData()
     NumberOfPoints = splinePoints.GetNumberOfPoints()
     for i in range(NumberOfPoints):
       p = splinePoints.GetPoint(i)
       pProj = [0, 0, 0]
       plane.ProjectPoint(p, pProj)
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
     TransformedPlane.SetInputConnection(reversePlane.GetOutputPort())
     TransformedPlane.SetTransform(splineTransform)

     finalModel = modelsLogic.AddModel(TransformedPlane.GetOutputPort())
     finalModel.GetDisplayNode().SetVisibility(False)
     finalModel.SetName("transformedPlane")

     implictSplinePlane = vtk.vtkImplicitPolyDataDistance()
     implictSplinePlane.SetInput(TransformedPlane.GetOutput())


     ###NOTE: clipping with the selection loop does not work as it will create an open surface... try to clip with the spline
     loop = vtk.vtkImplicitSelectionLoop()
     loop.SetLoop(PointsPolyData.GetPoints())
     v1 = -1 * plane.GetNormal()[0]
     v2 = -1 * plane.GetNormal()[1]
     v3 = -1 * plane.GetNormal()[2]
     loop.SetNormal(plane.GetNormal())

     # Clip the clipped input model with the spline plane so that only
     # parts of the scan about the spline plane are kept
     clippedInputWithLoop = vtk.vtkClipPolyData()
     clippedInputWithLoop.SetClipFunction(implictSplinePlane)  # should be loop
     clippedInputWithLoop.SetInputData(InputModel)
     clippedInputWithLoop.SetInsideOut(False)
     clippedInputWithLoop.Update()

     # No use the vtkPolyDataConnectivityFilter to extract the largest region
     connectedInput = vtk.vtkPolyDataConnectivityFilter()
     connectedInput.SetInputConnection(clippedInputWithLoop.GetOutputPort())
     connectedInput.SetExtractionModeToLargestRegion()
     connectedInput.Update()


     finalModel = modelsLogic.AddModel(connectedInput.GetOutputPort())
     finalModel.GetDisplayNode().SetVisibility(True)
     finalModel.SetName("clippedInput")
     finalModel.GetDisplayNode().BackfaceCullingOff()


     extrudeInputWithLoop = vtk.vtkLinearExtrusionFilter()
     extrudeInputWithLoop.SetInputConnection(connectedInput.GetOutputPort())
     extrudeInputWithLoop.SetScaleFactor(100)
     extrudeInputWithLoop.CappingOn()
     normVec = plane.GetNormal()
     if LeftBreast == True:
       extrudeInputWithLoop.SetVector((normVec[0]), (normVec[1]), (normVec[2]))
     else:
       extrudeInputWithLoop.SetVector((normVec[0] * -1), (normVec[1] * -1), (normVec[2] * -1))
     extrudeInputWithLoop.Update()

     # Do not use autoorient normals here, it will cause some surfaces not to be closed when
     # clipped with plane
     extrudeNormals = vtk.vtkPolyDataNormals()
     extrudeNormals.SetInputConnection(extrudeInputWithLoop.GetOutputPort())
     extrudeNormals.ComputePointNormalsOn()
     extrudeNormals.Update()

     clipped = vtk.vtkClipPolyData()
     clipped.SetClipFunction(loop)  # should be loop
     clipped.SetInputConnection(TransformedPlane.GetOutputPort())
     clipped.SetInsideOut(True)
     clipped.Update()

     clippedInputhWithPlane = vtk.vtkClipPolyData()
     clippedInputhWithPlane.SetClipFunction(implictSplinePlane)
     clippedInputhWithPlane.SetInputConnection(extrudeInputWithLoop.GetOutputPort())
     #When the function is not clipping correctly ie clipping behind breast instead
     # of front set inside out to False
     clippedInputhWithPlane.SetInsideOut(False)
     clippedInputhWithPlane.Update()

     clippedNormals = vtk.vtkPolyDataNormals()
     clippedNormals.SetInputConnection(clipped.GetOutputPort())
     clippedNormals.ComputePointNormalsOn()
     clippedNormals.ConsistencyOn()
     #If line 509 is changed to true then line 517 must also be commented out
     clippedNormals.FlipNormalsOn()
     clippedNormals.Update()

     implictInput = vtk.vtkImplicitPolyDataDistance()
     implictInput.SetInput(extrudeNormals.GetOutput())

     clippedWithImplictInput = vtk.vtkClipPolyData()
     clippedWithImplictInput.SetClipFunction(implictInput)  # should be loop
     clippedWithImplictInput.SetInputConnection(TransformedPlane.GetOutputPort())
     clippedWithImplictInput.SetInsideOut(True)
     clippedWithImplictInput.Update()

     appendClosedBreast = vtk.vtkAppendPolyData()
     appendClosedBreast.AddInputData(clippedInputhWithPlane.GetOutput()) #need to change these variables
     appendClosedBreast.AddInputData(clippedWithImplictInput.GetOutput())
     appendClosedBreast.Update()

     connectedOutput = vtk.vtkPolyDataConnectivityFilter()
     connectedOutput.SetInputConnection(appendClosedBreast.GetOutputPort())
     connectedOutput.SetExtractionModeToLargestRegion()
     connectedOutput.Update()

     cleanClosedBreast = vtk.vtkCleanPolyData()
     cleanClosedBreast.SetInputData(appendClosedBreast.GetOutput())
     cleanClosedBreast.Update()

     # output.SetPolyDataConnection(appendClosedBreast.GetOutputPort())
     # output.GetModelDisplayNode().VisibilityOn()

     finalModel = modelsLogic.AddModel(cleanClosedBreast.GetOutput())
     finalModel.GetDisplayNode().SetVisibility(True)
     finalModel.SetName("Closed Breast")
     finalModel.GetDisplayNode().BackfaceCullingOff()

     print("New Volume and SurfaceA")
     massProperties = vtk.vtkMassProperties()
     massProperties.SetInputData(cleanClosedBreast.GetOutput())
     volumeMP = massProperties.GetVolume()
     volumeMP = volumeMP / 1000
     volumeMP = round(volumeMP, 2)
     volume.setText(volumeMP)
     print('Volume in cc')
     print(volumeMP)
     surfaceAreaMP = massProperties.GetSurfaceArea()
     surfaceAreaMP = surfaceAreaMP / 100
     surfaceAreaMP = round(surfaceAreaMP, 2)
     surfaceArea.setText(surfaceAreaMP)
     print('Surface Area in cm^2')
     print(surfaceAreaMP)

     # Ensure the final model is closed so the volume computation is correct
     featureEdges = vtk.vtkFeatureEdges()
     featureEdges.FeatureEdgesOff()
     featureEdges.BoundaryEdgesOn()
     featureEdges.NonManifoldEdgesOn()
     featureEdges.SetInputData(cleanClosedBreast.GetOutput())
     featureEdges.Update()
     numberOfOpenEdges = featureEdges.GetOutput().GetNumberOfCells()

     if (numberOfOpenEdges > 0):
       print("Surface is not closed")
       print(numberOfOpenEdges)

     else:
       print("surface is closed")

  def AddVolumeNode(self):
    #Create volume for scene 
    imageSize = [512, 512, 512]
    imageSpacing = [1.0, 1.0, 1.0]
    voxelType = vtk.VTK_UNSIGNED_CHAR
    # Create an empty image volume
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType, 1)
    thresholder = vtk.vtkImageThreshold()
    thresholder.SetInputData(imageData)
    thresholder.SetInValue(0)
    thresholder.SetOutValue(0)
    # Create volume node
    volumeNode = slicer.vtkMRMLScalarVolumeNode()
    volumeNode.SetSpacing(imageSpacing)
    volumeNode.SetImageDataConnection(thresholder.GetOutputPort())
    # Add volume to scene
    slicer.mrmlScene.AddNode(volumeNode)
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    colorNode = slicer.util.getNode('Grey')
    displayNode.SetAndObserveColorNodeID(colorNode.GetID())
    volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
    volumeNode.CreateDefaultStorageNode()
    volumeNode.SetName("Volume Node")

  def run(self, inputModel, fidList, LeftBreast, volume, surfaceArea, output, planeType):
    """
    Run the actual algorithm
    """
    #If flat cut plane is chosen
    if planeType == True:
        self.cropWithPlane(inputModel, fidList, LeftBreast, volume, surfaceArea, output)
        self.AddVolumeNode()
    #else curved cut plane is chosen
    else:
        self.cropWithCurve(inputModel, fidList, LeftBreast, volume, surfaceArea, output)
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