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
    This module takes a 3D surface scans of a woman's torso and and lists of fiducials which mark the
    boundaries of the breasts. The breasts are separated from the rest of the scan and the volume is computed.
    """
        self.parent.acknowledgementText = """
    This module was developed in the Perklab at Queen's University.
""" 

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

        self.optionalInputSelector = slicer.qMRMLNodeComboBox()
        self.optionalInputSelector.nodeTypes = ["vtkMRMLModelNode"]
        self.optionalInputSelector.selectNodeUponCreation = True
        self.optionalInputSelector.addEnabled = False
        self.optionalInputSelector.removeEnabled = False
        self.optionalInputSelector.noneEnabled = True
        self.optionalInputSelector.showHidden = False
        self.optionalInputSelector.showChildNodeTypes = False
        self.optionalInputSelector.setMRMLScene(slicer.mrmlScene)
        self.optionalInputSelector.setToolTip("Pick the input to the algorithm.")
        parametersFormLayout.addRow("Input Initial (Optional): ", self.optionalInputSelector)

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

        self.applyButton = qt.QPushButton("Apply")
        self.applyButton.enabled = False
        parametersFormLayout.addRow("Breast Volume Computations", self.applyButton)

        # Volume and distance Labels
        self.VolumeLabelLeft = qt.QLabel()
        self.VolumeLabelLeft.setText("Volume in cc")
        parametersFormLayout.addRow("Left Breast Volume: ", self.VolumeLabelLeft)

        self.VolumeLabelRight = qt.QLabel()
        self.VolumeLabelRight.setText("Volume in cc")
        parametersFormLayout.addRow("Right Breast Volume: ", self.VolumeLabelRight)

        self.VolumeLabelDifference = qt.QLabel()
        self.VolumeLabelDifference.setText("Volume Difference in cc")
        parametersFormLayout.addRow("Breast Volume Difference: ", self.VolumeLabelDifference)

        self.MeanDistanceLabel = qt.QLabel()
        self.MeanDistanceLabel.setText("Mean distance in mm")
        parametersFormLayout.addRow("Mean Distance After Registration: ", self.MeanDistanceLabel)

        # Checkboxes for different module input
        self.noBreastCheck = qt.QCheckBox("")
        self.noBreastCheck.setCheckState(False)
        parametersFormLayout.addRow("Check to create torso mesh (No breast): ", self.noBreastCheck)

        self.infoLabel = qt.QLabel()
        self.infoLabel.setText("")

        parametersFormLayout.addRow("Please see documentation for examples of the below check box use", self.infoLabel)

        self.reverseNormalCheckL = qt.QCheckBox("")
        self.reverseNormalCheckL.setCheckState(False)
        parametersFormLayout.addRow("Check to reverse normal of left breast: ", self.reverseNormalCheckL)

        self.reverseNormalCheckR = qt.QCheckBox("")
        self.reverseNormalCheckR.setCheckState(False)
        parametersFormLayout.addRow("Check to reverse normal of right breast: ", self.reverseNormalCheckR)

        self.setInSideOutCheckL = qt.QCheckBox("")
        self.setInSideOutCheckL.setCheckState(False)
        parametersFormLayout.addRow("Check to set clip function to inside out of left breast: ", self.setInSideOutCheckL)

        self.setInSideOutCheckR = qt.QCheckBox("")
        self.setInSideOutCheckR.setCheckState(False)
        parametersFormLayout.addRow("Check to set clip function to inside out of right breast: ", self.setInSideOutCheckR)

        #    connections
        self.optionalInputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputLFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
        self.inputRFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
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

        if self.inputLFiducialSelector.currentNode() == None and self.inputRFiducialSelector.currentNode() == None:
            logging.error('Please input fiducials')

        # If the left markups are missing
        elif self.inputLFiducialSelector.currentNode() == None:
            logging.error('Please input left fiducials')

            self.reverseNormal = False
            if self.reverseNormalCheckR.isChecked():
                self.reverseNormal = True

            self.setInsideOut = False
            if self.setInSideOutCheckR.isChecked():
                self.setInsideOut = True

            self.noBreast = False
            if self.noBreastCheck.isChecked():
                self.noBreast = True

            logic = BreastReconstructionLogic()
            volumeL = 0

            [volumeR, meanDis] = logic.run(self.optionalInputSelector.currentNode(), self.inputModelSelector.currentNode(), self.inputRFiducialSelector.currentNode(), False,
                                           self.VolumeLabelRight, None, self.reverseNormal, self.setInsideOut, self.noBreast, True)
            if self.optionalInputSelector.currentNode() != None:
                self.MeanDistanceLabel.setText(meanDis)
            else:
                self.MeanDistanceLabel.setText("Mean distance in mm")

            self.VolumeLabelLeft.setText(0)
            volumeDif = volumeL - volumeR
            self.VolumeLabelDifference.setText(volumeDif)

        # If the right markups are missing
        elif self.inputRFiducialSelector.currentNode() == None:
            logging.error('Please input right fiducials')

            self.reverseNormal = False
            if self.reverseNormalCheckL.isChecked():
                self.reverseNormal = True

            self.setInsideOut = False
            if self.setInSideOutCheckL.isChecked():
                self.setInsideOut = True

            self.noBreast = False
            if self.noBreastCheck.isChecked():
                self.noBreast = True

            logic = BreastReconstructionLogic()
            volumeR = 0

            [volumeL, meanDis] = logic.run(self.optionalInputSelector.currentNode(), self.inputModelSelector.currentNode(), self.inputLFiducialSelector.currentNode(), True,
                                           self.VolumeLabelLeft, None, self.reverseNormal, self.setInsideOut, self.noBreast, True)
            if self.optionalInputSelector.currentNode() != None:
                self.MeanDistanceLabel.setText(meanDis)
            else:
                self.MeanDistanceLabel.setText("Mean distance in mm")


            self.VolumeLabelRight.setText(0)
            volumeDif = volumeL - volumeR
            self.VolumeLabelDifference.setText(volumeDif)

        # If both markups are present
        else:
            self.noBreast = False
            if self.noBreastCheck.isChecked():
                self.noBreast = True

            logic = BreastReconstructionLogic()

            self.reverseNormal = False
            if self.reverseNormalCheckL.isChecked():
                self.reverseNormal = True

            self.setInsideOut = False
            if self.setInSideOutCheckL.isChecked():
                self.setInsideOut = True

            [volumeL, meanDisL] = logic.run(self.optionalInputSelector.currentNode(), self.inputModelSelector.currentNode(), self.inputLFiducialSelector.currentNode(), True,
                                            self.VolumeLabelLeft, None, self.reverseNormal, self.setInsideOut, self.noBreast, True)
            self.reverseNormal = False
            if self.reverseNormalCheckR.isChecked():
                self.reverseNormal = True

            self.setInsideOut = False
            if self.setInSideOutCheckR.isChecked():
                self.setInsideOut = True

            [volumeR, meanDisR] = logic.run(self.optionalInputSelector.currentNode(), self.inputModelSelector.currentNode(), self.inputRFiducialSelector.currentNode(), False,
                                            self.VolumeLabelRight, None, self.reverseNormal, self.setInsideOut, self.noBreast, False)
            # Checks if the models were registered
            # If so display the mean registration error between the two healthy breast
            # this value is displayed instead of the value between the post BCT breasts so the
            # health breasts should remain constant in size (assuming healthy breast is larger)
            if self.optionalInputSelector.currentNode() != None:
                if volumeL > volumeR:
                    self.MeanDistanceLabel.setText(meanDisL)
                else:
                    self.MeanDistanceLabel.setText(meanDisR)
            else:
                self.MeanDistanceLabel.setText("Mean distance in mm")
            volumeDif = volumeL - volumeR
            self.VolumeLabelDifference.setText(volumeDif)

class BreastReconstructionLogic(ScriptedLoadableModuleLogic):

    def FiducialsToPolyData(self, fiducials, polyData):
        #create polydata from fiducial list
        points = vtk.vtkPoints()
        n = fiducials.GetNumberOfFiducials()
        for fiducialIndex in range( 0, n ):
            p = [0, 0, 0]
            fiducials.GetNthFiducialPosition(fiducialIndex, p)
            points.InsertNextPoint(p)
        tempPolyData = vtk.vtkPolyData()
        tempPolyData.SetPoints(points)

        vertex = vtk.vtkVertexGlyphFilter()
        vertex.SetInputData(tempPolyData)
        vertex.Update()

        polyData.ShallowCopy(vertex.GetOutput())

    def LeastSquaresPlane(self, modelNode, PointsPolyData, plane):
        '''
        Create a plane of best fit for the data points
        using the least squares approach
        '''
        numberofpoints = PointsPolyData.GetNumberOfPoints()
        #Creating Numpy array of all points
        allpoints = np.zeros((numberofpoints, 3))
        for i in range(numberofpoints):
            allpoints[i][0] = PointsPolyData.GetPoint(i)[0]
            allpoints[i][1] = PointsPolyData.GetPoint(i)[1]
            allpoints[i][2] = PointsPolyData.GetPoint(i)[2]

        # Computing the center of the points
        center = np.mean(allpoints, axis=0)

        #Create parameters for least  squares solution
        A = np.ones((numberofpoints,3))
        A[:,0] = allpoints[:,0]
        A[:,1] = allpoints[:,1]
        B = np.zeros((numberofpoints,1))
        B[:,0] = allpoints[:,2]
        #Solve least squares equation
        (a, b, c), resid, rank, s = np.linalg.lstsq(A, B,rcond=-1)
        #Construct Normal of Plane
        normal = (a, b, -1)
        #Normalized the plane normal
        nn = np.linalg.norm(normal)
        normal = normal / nn

        #Set original and normal of plane
        plane.SetOrigin(center)
        plane.SetNormal(normal)

    def surfaceRegistration(self, initialModel, inputModel):
        # Registers the input scan to the initial scan where the breast boundaries were defined

        # Uses ICP to create the transform
        icpTransform = vtk.vtkIterativeClosestPointTransform()
        icpTransform.SetSource(initialModel.GetPolyData())
        icpTransform.SetTarget(inputModel.GetPolyData())
        icpTransform.GetLandmarkTransform().SetModeToRigidBody()
        icpTransform.SetMaximumNumberOfIterations(100)
        icpTransform.Inverse()
        icpTransform.Modified()
        icpTransform.Update()

        # Applies the transform
        transformNode = slicer.vtkMRMLLinearTransformNode()
        transformNode.SetAndObserveMatrixTransformToParent(icpTransform.GetMatrix())
        slicer.mrmlScene.AddNode(transformNode)

        inputModel.SetAndObserveTransformNodeID(transformNode.GetID())
        inputModel.HardenTransform()

        return transformNode

    def distanceAfterRegistration(self, initialModel, inputModel, fidList, transformNode):
        # Computes the distance between the models following registration
        # computes the distance between the entire models but if the code below is
        # uncommented then only the distance between the breasts will be computed
        breastBoundPolyData = vtk.vtkPolyData()
        self.FiducialsToPolyData(fidList, breastBoundPolyData)

        # Create plane of best fit from input breast boundary fiducials
        plane = vtk.vtkPlane()
        self.LeastSquaresPlane(inputModel, breastBoundPolyData, plane)

        # UNCOMMENT out lines between if you would like to compute distance between only the breasts and not the
        # the entire models
        # creates loop from the breast boundary points
        # breastBound = vtk.vtkImplicitSelectionLoop()
        # breastBound.SetLoop(breastBoundPolyData.GetPoints())
        #
        # iM = vtk.vtkClipPolyData()
        # iM.SetClipFunction(breastBound)
        # iM.SetInputData(initialModel.GetPolyData())
        # iM.SetInsideOut(False)
        # iM.Update()
        #
        # M = vtk.vtkClipPolyData()
        # M.SetClipFunction(breastBound)
        # M.SetInputData(modelNode.GetPolyData())
        # M.SetInsideOut(False)
        # M.Update()
        #
        # sourcePolyData = iM.GetOutput()
        # targetPolyData = M.GetOutput()
        # STOP commenting now

        #Compute the mean distance after registration
        # If above is uncommented comment out the two lines following
        sourcePolyData = initialModel.GetPolyData()
        targetPolyData = inputModel.GetPolyData()

        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        dist2 = vtk.mutable(0.0)
        locator = vtk.vtkCellLocator()
        locator.SetDataSet(targetPolyData)
        locator.SetNumberOfCellsPerBucket(1)
        locator.BuildLocator()

        totalDistance = 0.0

        sourcePoints = sourcePolyData.GetPoints()
        n = sourcePoints.GetNumberOfPoints()
        m = vtk.vtkMath()
        for sourcePointIndex in xrange(n):
            sourcePointPos = [0, 0, 0]
            sourcePoints.GetPoint(sourcePointIndex, sourcePointPos)
            transformedSourcePointPos = [0, 0, 0, 1]
            sourcePointPos.append(1)
            transformNode.GetTransformToParent().MultiplyPoint(sourcePointPos, transformedSourcePointPos)
            surfacePoint = [0, 0, 0]
            transformedSourcePointPos.pop()
            locator.FindClosestPoint(transformedSourcePointPos, surfacePoint, cellId, subId, dist2)
            totalDistance = totalDistance + math.sqrt(dist2)

        return (round((totalDistance/n),2))

    def createTorsoModel(self, InputModel, loop):
        # creates a model node of the original scan where the breast has been removed
        modelsLogic = slicer.modules.models.logic()

        torso = vtk.vtkClipPolyData()
        torso.SetClipFunction(loop)
        torso.SetInputData(InputModel)
        torso.SetInsideOut(False)
        torso.Update()

        torsoModel = modelsLogic.AddModel(torso.GetOutput())
        torsoModel.GetDisplayNode().SetVisibility(True)
        torsoModel.SetName("noBreast")
        torsoModel.GetDisplayNode().BackfaceCullingOff()

    def createPlaneModel(self, InputModel, plane, breastFlag):
        #this function creates a model (visual representation) of the defined plane

        #The input is linearly extruded to create a closed input model so that when the cutter extracts the
        # region it is large enough to create a plane to cover the entire breast
        closedInputModel = vtk.vtkLinearExtrusionFilter()
        closedInputModel.SetInputData(InputModel)
        closedInputModel.SetScaleFactor(100)
        closedInputModel.CappingOn()
        closedInputModel.Update()

        clippedInput = vtk.vtkClipPolyData()
        clippedInput.SetInputData(closedInputModel.GetOutput())
        clippedInput.SetClipFunction(plane)
        clippedInput.SetValue(0)
        clippedInput.SetInsideOut(breastFlag)
        clippedInput.Update()

        cutterPlane = vtk.vtkCutter()
        cutterPlane.SetCutFunction(plane)
        cutterPlane.SetInputData(clippedInput.GetOutput())
        cutterPlane.Update()

        cutterModel = vtk.vtkPolyData()
        cutterModel = cutterPlane.GetOutput()
        surfPlane = vtk.vtkSurfaceReconstructionFilter()
        surfPlane.SetInputData(cutterModel)
        surfPlane.SetSampleSpacing(2.5)

        cfPlane = vtk.vtkContourFilter()
        cfPlane.SetInputConnection(surfPlane.GetOutputPort())
        cfPlane.SetValue(0, 0.0)
        reversePlane = vtk.vtkReverseSense()
        reversePlane.SetInputConnection(cfPlane.GetOutputPort())
        reversePlane.ReverseCellsOn()
        reversePlane.ReverseNormalsOn()

        return reversePlane

    def thinPlateSplineTransform(self, breastBoundPolyData,plane, modelNode):
        # Creates thin plate spline transform between target points
        # which are the breast boundary points (or code can be uncommented
        # to include more points) and the source points which are on
        # the plane of best fit

        # Target points using breast boundary points
        targetPoints = vtk.vtkPoints()
        for i in range(breastBoundPolyData.GetNumberOfPoints()):
            targetPoints.InsertNextPoint(breastBoundPolyData.GetPoint(i))

        # UNCOMMENT after this line to add more points

        # loop = vtk.vtkImplicitSelectionLoop()
        # loop.SetLoop(breastBoundPolyData.GetPoints())
        # torso = vtk.vtkClipPolyData()
        # torso.SetClipFunction(loop)
        # torso.SetInputData(modelNode.GetPolyData())
        # torso.SetInsideOut(False)
        # torso.Update()
        # o = plane.GetOrigin()
        # radius = 100
        # for i in range(torso.GetOutput().GetPoints().GetNumberOfPoints()): #can increment by larger number (ie 10) to speed up
        #     p = torso.GetOutput().GetPoints().GetPoint(i)
        #     dis = vtk.vtkMath.Distance2BetweenPoints(p,o)
        #     dis = math.sqrt(dis)
        #     if dis < radius: #This can be changed to input more or less of the input model
        #         targetPoints.InsertNextPoint(p)

        # END of section to add more points

        # Create source point set for thin plate spline by getting points on plane
        sourcePoints = vtk.vtkPoints()
        sourcePointsPolyData = vtk.vtkPolyData()
        NumberOfPoints = targetPoints.GetNumberOfPoints()
        for i in range(NumberOfPoints):
            p = targetPoints.GetPoint(i)
            pProj = [0, 0, 0]
            plane.ProjectPoint(p, pProj)
            sourcePoints.InsertNextPoint(pProj)
        sourcePointsPolyData.SetPoints(sourcePoints)

        # creates a thinPlate spline transform between the target and source points
        splineTransform = vtk.vtkThinPlateSplineTransform()
        splineTransform.SetSourceLandmarks(sourcePoints)
        splineTransform.SetTargetLandmarks(targetPoints)
        splineTransform.SetBasisToR()  # Since our points are 3D

        return splineTransform

    def createBreastModel(self, breastBound, InputModel, breastFlag, reverseNormal, setInsideOut, plane, planeWithSplineTransform):
        # Creates the cropped breast model from the original input model and the created warped posterior wall of the breast

        # Clip the input model with the breast boundary loop to isolate the breast
        clippedBreastModel = vtk.vtkClipPolyData()
        clippedBreastModel.SetClipFunction(breastBound)
        clippedBreastModel.SetInputData(InputModel)
        # This value below may need to be changed
        clippedBreastModel.SetInsideOut(True)
        clippedBreastModel.Update()

        # Now use the vtkPolyDataConnectivityFilter to extract the largest region
        connectedBreastModel = vtk.vtkPolyDataConnectivityFilter()
        connectedBreastModel.SetInputConnection(clippedBreastModel.GetOutputPort())
        connectedBreastModel.SetExtractionModeToLargestRegion()
        connectedBreastModel.Update()

        # linearly extrude the breast model to ensure the final surface will be airtight
        extrudeBreastModel = vtk.vtkLinearExtrusionFilter()
        extrudeBreastModel.SetInputConnection(connectedBreastModel.GetOutputPort())
        extrudeBreastModel.SetScaleFactor(100)
        extrudeBreastModel.CappingOn()
        normVec = plane.GetNormal()

        if breastFlag == True:
            normVec = [normVec[0] * -1, normVec[1] * -1, normVec[2] * -1]
        if reverseNormal == True:
            normVec = [normVec[0] * -1, normVec[1] * -1, normVec[2] * -1]
        extrudeBreastModel.SetVector(normVec)
        extrudeBreastModel.Update()

        # create implicit representation of the plane transformed with the spline
        implictSplinePlane = vtk.vtkImplicitPolyDataDistance()
        implictSplinePlane.SetInput(planeWithSplineTransform.GetOutput())

        # Now re-clip the now air-tight breast model wit the curved surface
        breastModel = vtk.vtkClipPolyData()
        breastModel.SetClipFunction(implictSplinePlane)
        breastModel.SetInputConnection(extrudeBreastModel.GetOutputPort())
        # This value below may need to be changed
        if setInsideOut == True:
            breastModel.SetInsideOut(True)
        else:
            breastModel.SetInsideOut(False)
        breastModel.Update()


        return breastModel.GetOutput()

    def createClippedPlane(self,breastBound, planeWithSplineTransform):
        # Takes the warped plane and crops it to become posterior wall of the breast

        clippedPlaneWithTransform = vtk.vtkClipPolyData()
        clippedPlaneWithTransform.SetClipFunction(breastBound)
        clippedPlaneWithTransform.SetInputConnection(planeWithSplineTransform.GetOutputPort())
        # This value below may need to be changed
        clippedPlaneWithTransform.SetInsideOut(True)
        clippedPlaneWithTransform.Update()

        connectedClippedPlaneWithTransform = vtk.vtkPolyDataConnectivityFilter()
        connectedClippedPlaneWithTransform.SetInputConnection(clippedPlaneWithTransform.GetOutputPort())
        connectedClippedPlaneWithTransform.SetExtractionModeToLargestRegion()
        connectedClippedPlaneWithTransform.Update()

        return connectedClippedPlaneWithTransform.GetOutput()

    def cropWithCurve(self, modelNode, fidList, breastFlag, reverseNormal, setInsideOut, torsoFlag):
        # This method takes in a surface scan and list of fiducials as input and computes the breast volume
        # A curved posterior chest wall is constructed to create the closed breast

        modelsLogic = slicer.modules.models.logic()
        #Check which breast volume is being computed for
        if breastFlag == True:
            name = "ClosedLeftBreast"
        else:
            name = "ClosedRightBreast"

        # Define parameters

        InputModel = modelNode.GetPolyData()
        breastBoundPolyData = vtk.vtkPolyData()
        self.FiducialsToPolyData(fidList, breastBoundPolyData)

        # Create plane of best fit from input breast boundary fiducials
        plane = vtk.vtkPlane()
        self.LeastSquaresPlane(modelNode, breastBoundPolyData, plane)

        #creates loop from the breast boundary points
        breastBound = vtk.vtkImplicitSelectionLoop()
        breastBound.SetLoop(breastBoundPolyData.GetPoints())

        #creates the torso model when the torso flag is set
        if torsoFlag == True:
            self.createTorsoModel(InputModel,breastBound)

        # creates model for the plane of best fit
        planeModel = vtk.vtkReverseSense()
        planeModel = self.createPlaneModel(InputModel, plane, breastFlag)

        splineTransform = vtk.vtkThinPlateSplineTransform()
        splineTransform = self.thinPlateSplineTransform(breastBoundPolyData,plane, modelNode)

        #apply spline transform to the plane visualization
        planeWithSplineTransform = vtk.vtkTransformPolyDataFilter()
        planeWithSplineTransform.SetInputConnection(planeModel.GetOutputPort())
        planeWithSplineTransform.SetTransform(splineTransform)

        #create model of the transformed plane
        transformedPlaneModel = modelsLogic.AddModel(planeWithSplineTransform.GetOutputPort())
        transformedPlaneModel.GetDisplayNode().SetVisibility(False)
        transformedPlaneModel.SetName("transformedPlane")

        # Creates cropped breast model
        breastModel = self.createBreastModel(breastBound, InputModel, breastFlag, reverseNormal, setInsideOut, plane, planeWithSplineTransform)

        # Creates cropped-posterior wall
        posteriorWallModel = vtk.vtkPolyData()
        posteriorWallModel = self.createClippedPlane(breastBound,planeWithSplineTransform)

        # Creates closed breast model
        appendClosedBreast = vtk.vtkAppendPolyData()
        appendClosedBreast.AddInputData(breastModel) #Breasts
        appendClosedBreast.AddInputData(posteriorWallModel) #Chest Wall
        appendClosedBreast.Update()

        # Applies cleaning filter to closed breast model
        cleanClosedBreast = vtk.vtkCleanPolyData()
        cleanClosedBreast.SetInputData(appendClosedBreast.GetOutput())
        cleanClosedBreast.Update()

        # Added closed breast model to slicer models
        closedBreastModel = modelsLogic.AddModel(cleanClosedBreast.GetOutput())
        closedBreastModel.GetDisplayNode().SetVisibility(True)
        closedBreastModel.SetName(name)
        closedBreastModel.GetDisplayNode().BackfaceCullingOff()

        slicer.mrmlScene.RemoveNode(transformedPlaneModel)

        return closedBreastModel

    def AddVolumeNode(self):
        # Create empty volume node
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
        #This function creates a segmentation from the output model (ie cropped breast models)
        # and computed the volume of the segmentation, this is done to ensure the volume is computed correctly

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


    def run(self, initialModel, inputModel, fidList, LeftBreast, volume, modelFids, reverseNormal,
            setInsideOut, noBreast, registrationFlag):
        """
        Run the actual algorithm
        """
        meanDis = 0
        transformNode = slicer.vtkMRMLLinearTransformNode()
        if initialModel != None:
        #If the models need to registered
            if registrationFlag == True:
                transformNode = self.surfaceRegistration(initialModel, inputModel)
                meanDis = self.distanceAfterRegistration(initialModel, inputModel, fidList, transformNode)
            else:
                meanDis = self.distanceAfterRegistration(initialModel, inputModel, fidList, transformNode)
        # Creates the closed breast model
        modelNode = self.cropWithCurve(inputModel, fidList, LeftBreast, reverseNormal, setInsideOut, noBreast)
        self.AddVolumeNode()
        # Computes the volume of the closed breast model
        volumeVal = self.AddSegmentNode(modelNode, volume)
        return [volumeVal, meanDis]

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

