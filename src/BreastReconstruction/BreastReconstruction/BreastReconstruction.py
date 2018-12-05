import os
import unittest
import vtk, qt, ctk, slicer
import math
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
from random import *
import vtkSegmentationCorePython as vtkSegmentationCore

#
# BreastReconstruction
#

class BreastReconstruction(ScriptedLoadableModule):

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "BreastReconstruction"
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

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

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

        #input left Fiducal selector
        self.inputLFiducialSelector = slicer.qSlicerSimpleMarkupsWidget()
        self.inputLFiducialSelector.tableWidget().hide()
        self.inputLFiducialSelector.setMRMLScene(slicer.mrmlScene)
        self.inputLFiducialSelector.setToolTip("Pick the fiducials to define the left breast region of interest.")
        self.inputLFiducialSelector.setNodeBaseName ("LeftFiducials")

        parametersFormLayout.addRow("Left fiducials: ", self.inputLFiducialSelector)
        # Enable place multiple markups by default
        placeWidget = self.inputLFiducialSelector.markupsPlaceWidget()
        placeWidget.placeMultipleMarkups = slicer.qSlicerMarkupsPlaceWidget.ForcePlaceMultipleMarkups
        placeWidget.placeModeEnabled = False
        placeWidget.placeModeEnabled = True

        # input right  Fiducal selector
        self.inputRFiducialSelector = slicer.qSlicerSimpleMarkupsWidget()
        self.inputRFiducialSelector.tableWidget().hide()
        self.inputRFiducialSelector.setMRMLScene(slicer.mrmlScene)
        self.inputRFiducialSelector.setToolTip("Pick the fiducials to define the right breast region of interest.")
        self.inputRFiducialSelector.setNodeBaseName("RightFiducials")

        parametersFormLayout.addRow("Right fiducials: ", self.inputRFiducialSelector)

        placeWidget = self.inputRFiducialSelector.markupsPlaceWidget()
        placeWidget.placeMultipleMarkups = slicer.qSlicerMarkupsPlaceWidget.ForcePlaceMultipleMarkups
        placeWidget.placeModeEnabled = False
        placeWidget.placeModeEnabled = True

        #Add buttons
        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.enabled = False
        parametersFormLayout.addRow("Breast Volume Computations", self.applyButton)

        self.VolumeLabelLeft = qt.QLabel()
        self.VolumeLabelLeft.setText("Volume in cc")
        parametersFormLayout.addRow("Left Breast Volume: ", self.VolumeLabelLeft)

        self.VolumeLabelRight = qt.QLabel()
        self.VolumeLabelRight.setText("Volume in cc")
        parametersFormLayout.addRow("Right Breast Volume: ", self.VolumeLabelRight)

        self.VolumeLabelDifference = qt.QLabel()
        self.VolumeLabelDifference.setText("Volume Difference in cc")
        parametersFormLayout.addRow("Breast Volume Difference: ", self.VolumeLabelDifference)

        self.error1Check = qt.QCheckBox("Error Type 1")
        self.error1Check.setCheckState(False)
        parametersFormLayout.addWidget(self.error1Check)

        self.error2Check = qt.QCheckBox("Error Type 2")
        self.error2Check.setCheckState(False)
        parametersFormLayout.addWidget(self.error2Check)

        #    connections
        self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputLFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputRFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        #self.FiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.applyButton.connect('clicked(bool)', self.onApplyButton)

        # Add vertical spacer
        self.layout.addStretch(1)

        # Refresh Apply button state
        self.onSelect()

    def cleanup(self):
        pass

    def onSelect(self):
        self.applyButton.enabled = self.inputModelSelector.currentNode()

    def onApplyButton(self):

        if self.inputLFiducialSelector.currentNode() == None:
            logging.error('Please enter input left fiducials')

        if self.inputRFiducialSelector.currentNode() == None:
            logging.error('Please enter input right fiducials')
        elif self.inputModelSelector == None:
            self.error("Please enter input model")
        else:
            self.error1 = False
            self.error2 = False
            if self.error1Check.isChecked():
                self.error1 = True
            if self.error2Check.isChecked():
                self.error2 = True
            logic = BreastReconstructionLogic()
            volumeL = logic.run(self.inputModelSelector.currentNode(), self.inputLFiducialSelector.currentNode(), True,
                      self.VolumeLabelLeft, None, self.error1, self.error2)
            volumeR = logic.run(self.inputModelSelector.currentNode(), self.inputRFiducialSelector.currentNode(), False,
                      self.VolumeLabelRight, None, self.error1, self.error2)
            volumeDif = volumeL - volumeR
            self.VolumeLabelDifference.setText(volumeDif)
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

    def LeastSquaresPlane(self, modelNode, PointsPolyData, plane):
        '''
        create a least square fitted plane to the points
        create plane by finding center of mass and direction vector
        '''
        #Store fiducials as PolyData
        #PointsPolyData = vtk.vtkPolyData()
        #self.FiducialsToPolyData(fidList, PointsPolyData)
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

    #This is the method used to compute breast volume
    def cropWithCurve(self, modelNode, fidList, LeftBreast, volume, modelFids, error1, error2):
        modelsLogic = slicer.modules.models.logic()
        #Check which breast volume is being computed for
        if LeftBreast == True:
            name = "ClosedLeftBreast"
        else:
            name = "ClosedRightBreast"

        PointsPolyData = vtk.vtkPolyData()
        #This line was used to introduce random error into the fiduical placement
        #self.errorSim(fidList, PointsPolyData)
        self.FiducialsToPolyData(fidList, PointsPolyData)

        # Clip the input model with the plane defined by input points
        plane = vtk.vtkPlane()
        self.LeastSquaresPlane(modelNode, PointsPolyData, plane)
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

        # #****************************************************************************
        # spline model to model back of the chest wall
        splinePoints = vtk.vtkPoints()

        #these commented out lines where for if the user inputed seperate points to select points of model surface istead of just the selection loop
        # for i in range(modelPoints.GetNumberOfPoints()):
        #     splinePoints.InsertNextPoint(modelPoints.GetPoint(i))
        for i in range(PointsPolyData.GetNumberOfPoints()):
            splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))
        #creates implicit representation of spline
        spline = vtk.vtkParametricSpline()
        spline.SetPoints(splinePoints)
        spline.ClosedOn()
        functionSource = vtk.vtkParametricFunctionSource()
        functionSource.SetParametricFunction(spline)
        implictSpline = vtk.vtkImplicitPolyDataDistance()
        implictSpline.SetInput(functionSource.GetOutput())

        # Create a new point set which is the contour points projected onto the plane
        projectedPoints = vtk.vtkPoints()
        projectedPointsPolyData = vtk.vtkPolyData()
        NumberOfPoints = splinePoints.GetNumberOfPoints()
        for i in range(NumberOfPoints):
            p = splinePoints.GetPoint(i)
            pProj = [0, 0, 0]
            plane.ProjectPoint(p, pProj)
            projectedPoints.InsertNextPoint(pProj)
        projectedPointsPolyData.SetPoints(projectedPoints)

        #creating a spline curve and implicit representation of projectedPoints
        spline2 = vtk.vtkParametricSpline()
        spline2.SetPoints(projectedPointsPolyData.GetPoints())
        spline2.ClosedOn()
        functionSource2 = vtk.vtkParametricFunctionSource()
        functionSource2.SetParametricFunction(spline2)

        #creates a thinspline trasnform between the original points and the projectedPoints
        splineTransform = vtk.vtkThinPlateSplineTransform()
        splineTransform.SetSourceLandmarks(projectedPoints)
        splineTransform.SetTargetLandmarks(PointsPolyData.GetPoints())
        splineTransform.SetBasisToR() #Since our points are 3D

        #apply spline transform to the plane visualization
        planeWithSplineTransform = vtk.vtkTransformPolyDataFilter()
        planeWithSplineTransform.SetInputConnection(reversePlane.GetOutputPort())
        planeWithSplineTransform.SetTransform(splineTransform)

        #create model of the transformed plane
        transformedPlaneModel = modelsLogic.AddModel(planeWithSplineTransform.GetOutputPort())
        transformedPlaneModel.GetDisplayNode().SetVisibility(False)
        transformedPlaneModel.SetName("transformedPlane")

        #create implicit representation of the plane transformed with the spline
        implictSplinePlane = vtk.vtkImplicitPolyDataDistance()
        implictSplinePlane.SetInput(planeWithSplineTransform.GetOutput())


        loop = vtk.vtkImplicitSelectionLoop()
        loop.SetLoop(PointsPolyData.GetPoints())
        v1 = -1 * plane.GetNormal()[0]
        v2 = -1 * plane.GetNormal()[1]
        v3 = -1 * plane.GetNormal()[2]
        loop.SetNormal(plane.GetNormal())

        clippedInputWithLoop = vtk.vtkClipPolyData()
        clippedInputWithLoop.SetClipFunction(loop)
        clippedInputWithLoop.SetInputData(InputModel)
        #####Change here as well################################
        clippedInputWithLoop.SetInsideOut(True)
        clippedInputWithLoop.Update()

        # No use the vtkPolyDataConnectivityFilter to extract the largest region
        connectedInput = vtk.vtkPolyDataConnectivityFilter()
        connectedInput.SetInputConnection(clippedInputWithLoop.GetOutputPort())
        connectedInput.SetExtractionModeToLargestRegion()
        connectedInput.Update()

        extrudeInputWithLoop = vtk.vtkLinearExtrusionFilter()
        extrudeInputWithLoop.SetInputConnection(connectedInput.GetOutputPort())
        extrudeInputWithLoop.SetScaleFactor(100)
        extrudeInputWithLoop.CappingOn()
        normVec = plane.GetNormal()
        if (LeftBreast == True) or (error2 == True):
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

        clippedTransformedPlane = vtk.vtkClipPolyData()
        clippedTransformedPlane.SetClipFunction(loop)  # should be loop
        clippedTransformedPlane.SetInputConnection(planeWithSplineTransform.GetOutputPort())
        clippedTransformedPlane.SetInsideOut(True)
        clippedTransformedPlane.Update()

        clippedInputhWithPlane = vtk.vtkClipPolyData()
        clippedInputhWithPlane.SetClipFunction(implictSplinePlane)
        clippedInputhWithPlane.SetInputConnection(extrudeInputWithLoop.GetOutputPort())
        #When the function is not clipping correctly ie clipping behind breast instead
        # of front set inside out to False
        ########Change this line below
        #************************************************************

        clippedInputhWithPlane.SetInsideOut(False)
        if error1 == True or error2 == True:
            clippedInputhWithPlane.SetInsideOut(True)
        clippedInputhWithPlane.Update()

        clippedTransformedPlaneNormals = vtk.vtkPolyDataNormals()
        clippedTransformedPlaneNormals.SetInputConnection(clippedTransformedPlane.GetOutputPort())
        clippedTransformedPlaneNormals.ComputePointNormalsOn()
        clippedTransformedPlaneNormals.ConsistencyOn()
        #If line 605 is changed to true then line 596 must also be commented out
        ####change this line below
        #*******************************************
        clippedTransformedPlaneNormals.FlipNormalsOn()
        clippedTransformedPlaneNormals.Update()

        implictInput = vtk.vtkImplicitPolyDataDistance()
        implictInput.SetInput(extrudeNormals.GetOutput())

        clippedPlaneWithTransform = vtk.vtkClipPolyData()
        clippedPlaneWithTransform.SetClipFunction(implictInput)
        clippedPlaneWithTransform.SetInputConnection(planeWithSplineTransform.GetOutputPort())
        ###this is line must also be changed###
        clippedPlaneWithTransform.SetInsideOut(True)
        clippedPlaneWithTransform.Update()

        connectedClippedPlaneWithTransform = vtk.vtkPolyDataConnectivityFilter()
        connectedClippedPlaneWithTransform.SetInputConnection(clippedPlaneWithTransform.GetOutputPort())
        connectedClippedPlaneWithTransform.SetExtractionModeToLargestRegion()
        connectedClippedPlaneWithTransform.Update()

        appendClosedBreast = vtk.vtkAppendPolyData()
        appendClosedBreast.AddInputData(clippedInputhWithPlane.GetOutput()) #need to change these variables
        appendClosedBreast.AddInputData(connectedClippedPlaneWithTransform.GetOutput())
        appendClosedBreast.Update()

        cleanClosedBreast = vtk.vtkCleanPolyData()
        cleanClosedBreast.SetInputData(appendClosedBreast.GetOutput())
        cleanClosedBreast.Update()

        closedBreastModel = modelsLogic.AddModel(cleanClosedBreast.GetOutput())
        closedBreastModel.GetDisplayNode().SetVisibility(True)
        closedBreastModel.SetName(name)
        closedBreastModel.GetDisplayNode().BackfaceCullingOff()

        slicer.mrmlScene.RemoveNode(transformedPlaneModel)

        return closedBreastModel

    def AddVolumeNode(self):
        # Create volume node
        volumeNode = slicer.mrmlScene.GetNthNodeByClass(1, 'vtkMRMLScalarVolumeNode')
        if volumeNode == None:
            volumeNode = slicer.vtkMRMLScalarVolumeNode()
            slicer.mrmlScene.AddNode(volumeNode)

        displayNode = slicer.mrmlScene.GetNthNodeByClass(1,'vtkMRMLScalarVolumeDisplayNode')
        if displayNode == None:
            displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
        colorNode = slicer.util.getNode('Grey')
        displayNode.SetAndObserveColorNodeID(colorNode.GetID())
        volumeNode.SetAndObserveDisplayNodeID(displayNode.GetID())
        volumeNode.CreateDefaultStorageNode()
        volumeNode.SetName("Volume Node")

    def AddSegmentNode(self, modelNode, volumeLabel):
        #This function creates a segmentation from the output model (ie cropped breast models) and computed the volume

        segmentationNode = slicer.vtkMRMLSegmentationNode()
        slicer.mrmlScene.AddNode(segmentationNode)

        displayNode = slicer.mrmlScene.GetNthNodeByClass(1, 'vtkMRMLScalarVolumeDisplayNode')
        if displayNode == None:
            displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()

        segmentationNode.SetAndObserveDisplayNodeID(displayNode.GetID())
        segmentationNode.SetName("Model Segmentation Node")

        slicer.modules.segmentations.logic().ImportModelToSegmentationNode(modelNode, segmentationNode)
        segment = segmentationNode.GetSegmentation().GetNthSegment(0)
        segBinaryLabelName = vtkSegmentationCore.vtkSegmentationConverter.GetSegmentationBinaryLabelmapRepresentationName()
        segmentLabelmap = segment.GetRepresentation(segBinaryLabelName)
        # We need to know exactly the value of the segment voxels, apply threshold to make force the selected label value
        labelValue = 1
        backgroundValue = 0
        thresh = vtk.vtkImageThreshold()
        thresh.SetInputData(segmentLabelmap)
        thresh.ThresholdByLower(0)
        thresh.SetInValue(backgroundValue)
        thresh.SetOutValue(labelValue)
        thresh.SetOutputScalarType(vtk.VTK_UNSIGNED_CHAR)
        thresh.Update()

        #  Use binary labelmap as a stencil
        stencil = vtk.vtkImageToImageStencil()
        stencil.SetInputData(thresh.GetOutput())
        stencil.ThresholdByUpper(labelValue)
        stencil.Update()

        stat = vtk.vtkImageAccumulate()
        stat.SetInputData(thresh.GetOutput())
        stat.SetStencilData(stencil.GetOutput())
        stat.Update()

        # Add data to statistics list
        cubicMMPerVoxel = reduce(lambda x, y: x * y, segmentLabelmap.GetSpacing())
        ccPerCubicMM = 0.001
        stats = {}
        volume = round(stat.GetVoxelCount() * cubicMMPerVoxel * ccPerCubicMM,2)
        volumeLabel.setText(volume)
        slicer.mrmlScene.RemoveNode(segmentationNode)
        return volume


    def run(self, inputModel, fidList, LeftBreast, volume, modelFids, error1, error2):
        """
        Run the actual algorithm
        """
        modelNode = self.cropWithCurve(inputModel, fidList, LeftBreast, volume, modelFids, error1, error2)
        self.AddVolumeNode()
        volumeVal = self.AddSegmentNode(modelNode, volume)
        return volumeVal

        logging.info('Processing completed')



class BreastReconstructionTest(ScriptedLoadableModuleTest):

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

        self.delayDisplay("Starting the test")
        self.delayDisplay('Test passed!')