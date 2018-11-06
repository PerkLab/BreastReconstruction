import os
import unittest
import vtk, qt, ctk, slicer
import math
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
from random import *

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

        # # points on Model
        # self.FiducialSelector = slicer.qSlicerSimpleMarkupsWidget()
        # self.FiducialSelector.tableWidget().hide()
        # self.FiducialSelector.setMRMLScene(slicer.mrmlScene)
        # self.FiducialSelector.setToolTip("Points on Model.")
        # self.FiducialSelector.setNodeBaseName("PointFiducials")
        #
        # parametersFormLayout.addRow("Point fiducials: ", self.FiducialSelector)
        # Enable place multiple markups by default
        placeWidget = self.inputRFiducialSelector.markupsPlaceWidget()
        placeWidget.placeMultipleMarkups = slicer.qSlicerMarkupsPlaceWidget.ForcePlaceMultipleMarkups
        placeWidget.placeModeEnabled = False
        placeWidget.placeModeEnabled = True

        self.planeCheckedBox = qt.QCheckBox("Plane cut")
        self.planeCheckedBox.setCheckState(False)
        parametersFormLayout.addWidget(self.planeCheckedBox)

        self.curvedCheckedBox = qt.QCheckBox("Curve cut")
        self.curvedCheckedBox.setCheckState(True)
        parametersFormLayout.addWidget(self.curvedCheckedBox)

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

        self.SurfaceAreaLabelLeft = qt.QLabel()
        self.SurfaceAreaLabelLeft.setText("Surface Area in cm^2")
        parametersFormLayout.addRow("Left Breast Surface Area: ", self.SurfaceAreaLabelLeft)

        self.SurfaceAreaLabelRight = qt.QLabel()
        self.SurfaceAreaLabelRight.setText("Surface Area in cm^2")
        parametersFormLayout.addRow("Right Breast Surface Area: ", self.SurfaceAreaLabelRight)

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
        if self.planeCheckedBox.isChecked():
            #Use flat cut plane
            planeType = True
        elif self.curvedCheckedBox.isChecked():
            #Use curved cut plane
            planeType = False
        else:
            logging.error('Please select either plane or curved cut')

        if self.inputLFiducialSelector.currentNode() == None:
            logging.error('Please enter input left fiducials')

        if self.inputRFiducialSelector.currentNode() == None:
            logging.error('Please enter input right fiducials')
        elif self.inputModelSelector == None:
            self.error("Please enter input model")
        else:
            logic = BreastReconstructionLogic()
            volumeL = logic.run(self.inputModelSelector.currentNode(), self.inputLFiducialSelector.currentNode(), True,
                      self.VolumeLabelLeft, self.SurfaceAreaLabelLeft,
                      planeType, None)#self.FiducialSelector.currentNode())
            volumeR = logic.run(self.inputModelSelector.currentNode(), self.inputRFiducialSelector.currentNode(), False,
                      self.VolumeLabelRight, self.SurfaceAreaLabelRight,
                      planeType, None)#self.FiducialSelector.currentNode())
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


    def cropWithPlane(self, modelNode, fidList, LeftBreast, volume, surfaceArea):

        # # Check which breast volume is being computed for
        # if LeftBreast == True:
        #     name = "ClosedLeftBreast"
        # else:
        #     name = "ClosedRightBreast"
        #
        # modelsLogic = slicer.modules.models.logic()
        # InputModel = modelNode.GetPolyData()
        #
        #
        # #Test creating parametric ellipsoid######################################
        # modelsLogic = slicer.modules.models.logic()
        # PointsPolyData = vtk.vtkPolyData()
        # self.FiducialsToPolyData(fidList, PointsPolyData)
        # NumberOfPoints = PointsPolyData.GetNumberOfPoints()
        # # Compute the center of mass of all points
        # CenterOfMass = vtk.vtkCenterOfMass()
        # CenterOfMass.SetInputData(PointsPolyData)
        # CenterOfMass.SetUseScalarsAsWeights(False)
        # CenterOfMass.Update()
        #
        # squaredDistance1 = vtk.vtkMath().Distance2BetweenPoints(CenterOfMass.GetCenter(), PointsPolyData.GetPoint(1))
        # dis1 = math.sqrt(squaredDistance1)
        #
        # squaredDistance2 = vtk.vtkMath().Distance2BetweenPoints(CenterOfMass.GetCenter(), PointsPolyData.GetPoint(2))
        # dis2 = math.sqrt(squaredDistance2)
        #
        # squaredDistance3 = vtk.vtkMath().Distance2BetweenPoints(CenterOfMass.GetCenter(), PointsPolyData.GetPoint(3))
        # dis3 = math.sqrt(squaredDistance3)
        #
        # squaredDistance4 = vtk.vtkMath().Distance2BetweenPoints(CenterOfMass.GetCenter(), PointsPolyData.GetPoint(4))
        # dis4 = math.sqrt(squaredDistance4)
        #
        # elipsoid = vtk.vtkParametricEllipsoid()
        # elipsoid.SetZRadius(0.5)
        # elipsoid.SetXRadius(max(dis1, dis3))
        # elipsoid.SetYRadius(max(dis2,dis4))
        #
        # sphere = vtk.vtkSphereSource()
        # sphere.SetThetaResolution(12)
        # sphere.SetPhiResolution(12)
        # sphere.SetCenter(CenterOfMass.GetCenter())
        # sphere.SetRadius(max(dis1,dis3))
        #
        # functionSource = vtk.vtkParametricFunctionSource()
        # functionSource.SetParametricFunction(elipsoid)
        # functionSource.Update()
        #
        # ##########################################################################################
        # eModel = modelsLogic.AddModel(functionSource.GetOutputPort())
        # eModel.GetDisplayNode().SetVisibility(True)
        # eModel.SetName("elipsoid")
        # eModel.GetDisplayNode().BackfaceCullingOff()
        #
        # # sModel = modelsLogic.AddModel(sphere.GetOutputPort())
        # # sModel.GetDisplayNode().SetVisibility(True)
        # # sModel.SetName("sphere")
        # # sModel.GetDisplayNode().BackfaceCullingOff()
        #
        #
        #
        #
        #
        # #####################################################################################################
        # #Creating best fit plane
        # plane = vtk.vtkPlane()
        # pointsPolyData = vtk.vtkPolyData()
        # self.FiducialsToPolyData(fidList, PointsPolyData)
        # NumberOfPoints = PointsPolyData.GetNumberOfPoints()
        # # Compute the center of mass of all points
        # CenterOfMass = vtk.vtkCenterOfMass()
        # CenterOfMass.SetInputData(PointsPolyData)
        # CenterOfMass.SetUseScalarsAsWeights(False)
        # CenterOfMass.Update()
        # center = CenterOfMass.GetCenter()
        #
        # tempPlane = vtk.vtkPlane()
        # tempPlane.SetOrigin(center)
        #
        # plane.SetOrigin(center)
        #
        # bestDistance = float('inf')
        # for i in range((NumberOfPoints)):
        #     for j in range(1,(NumberOfPoints)):
        #         for k in range(2, NumberOfPoints):
        #             if (i != j) and (i != k) and (j != k):
        #                 #Create plane of best fit from points
        #                 triangle = vtk.vtkTriangle()
        #                 pi = PointsPolyData.GetPoint(i)
        #                 pj = PointsPolyData.GetPoint(j)
        #                 pk = PointsPolyData.GetPoint(k)
        #                 # compute the normal of the 3 selected points
        #                 normalVector = [0.0, 0.0, 0.0]
        #                 triangle.ComputeNormal(pi, pj, pk, normalVector)
        #                 tempPlane.SetNormal(normalVector)
        #                 distance = 0
        #                 for p in range(NumberOfPoints):
        #                     # compute the distance from each point to the plane
        #                     point = PointsPolyData.GetPoint(p)
        #                     distance = distance + abs(tempPlane.DistanceToPlane(point))
        #                 averageDistance = distance / NumberOfPoints
        #                 if averageDistance < bestDistance:
        #                     # select plane if average distance value is lower
        #                     bestDistance = averageDistance
        #                     plane.SetNormal(normalVector)
        #
        # #create visual representation of the plane to add to the scene
        # cutterPlane = vtk.vtkCutter()
        # cutterPlane.SetCutFunction(plane)
        # cutterPlane.SetInputData(InputModel)
        # cutterPlane.Update()
        #
        # cutterModel = vtk.vtkPolyData()
        # cutterModel = cutterPlane.GetOutput()
        # surfPlane = vtk.vtkSurfaceReconstructionFilter()
        # surfPlane.SetInputData(cutterModel)
        # cfPlane = vtk.vtkContourFilter()
        # cfPlane.SetInputConnection(surfPlane.GetOutputPort())
        # cfPlane.SetValue(0, 0.0)
        # reversePlane = vtk.vtkReverseSense()
        # reversePlane.SetInputConnection(cfPlane.GetOutputPort())
        # reversePlane.ReverseCellsOn()
        # reversePlane.ReverseNormalsOn()
        #
        # pModel = modelsLogic.AddModel(reversePlane.GetOutputPort())
        # pModel.GetDisplayNode().SetVisibility(True)
        # pModel.SetName("plane")
        # pModel.GetDisplayNode().BackfaceCullingOff()
        #     #need to compute the transform between the plane and the ellipsoid
        #
        #
        # ##################################################################################
        # #Computing rotation between ellipsoid and plane by computing the rotation between the two normals
        # #Note: The ellipsoid is parallel to the x y axis, so it's normal vector is the z axis
        # mathvtk = vtk.vtkMath()
        # norm1 = [0, 0, 1]
        # norm2 = plane.GetNormal()
        # vec = [0,0,0]
        # mathvtk.Cross(norm2,norm1,vec)
        # cosT = mathvtk.Dot(norm2,norm1)
        # sinT = mathvtk.Norm(vec)
        # theta = math.atan2(sinT, cosT)
        # if sinT != 0:
        #     vec[0] /= sinT
        #     vec[1] /= sinT
        #     vec[2] /= sinT
        # #Convert to quaternion
        # cosT = math.cos(0.5*theta)
        # sinT = math.sin(0.5*theta)
        # quat = [0,0,0,0]
        # quat[0] = cosT
        # quat[1] = vec[0]*sinT
        # quat[2] = vec[1]*sinT
        # quat[3] = vec[2]*sinT
        # #convert to a matrix
        # mat = [[0,0,0],[0,0,0],[0,0,0]]
        # mathvtk.QuaternionToMatrix3x3(quat, mat)
        # cent = (CenterOfMass.GetCenter())
        # cent2 = [0, 0 ,0]
        # cent2[0] = cent[0]*-1
        # cent2[1] = cent[1] *-1
        # cent2[2] = cent[2] *-1
        # transMat = [mat[0][0], mat[1][0], mat[2][0],0,
        #             mat[0][1], mat[1][1], mat[2][1], 0,
        #             mat[0][2], mat[1][2], mat[2][2], 0,
        #             0,0,0,1]
        #             #cent[0], cent[1], cent[2], 1]
        #
        # logicTrans = slicer.vtkSlicerTransformLogic()
        # transformNode1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode")
        # eModel.SetAndObserveTransformNodeID(transformNode1.GetID())
        # transform = vtk.vtkTransform()
        # transform.SetMatrix(transMat)
        # #transform.Translate(cent2)
        # transformNode1.SetMatrixTransformToParent(transform.GetMatrix())
        #
        # logicTrans.hardenTransform(transformNode1)
        #
        # transformNode1 = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTransformNode")
        # eModel.SetAndObserveTransformNodeID(transformNode1.GetID())
        # transform = vtk.vtkTransform()
        # transform.Translate(cent)
        # transformNode1.SetMatrixTransformToParent(transform.GetMatrix())
        #
        # logicTrans.hardenTransform(transformNode1)


#############################################################

        # Check which breast volume is being computed for
        if LeftBreast == True:
          name = "ClosedLeftBreast"
        else:
          name = "ClosedRightBreast"

        modelsLogic = slicer.modules.models.logic()
        InputModel = modelNode.GetPolyData()

        PointsPolyData = vtk.vtkPolyData()
        #self.errorSim(fidList, PointsPolyData)
        self.FiducialsToPolyData(fidList, PointsPolyData)

        # Clip the input model with the plane defined by input points
        plane = vtk.vtkPlane()
        self.LeastSquaresPlane(modelNode, PointsPolyData, plane)


        # create visual representation of the plane to add to the scene
        cutterPlane = vtk.vtkCutter()
        cutterPlane.SetCutFunction(plane)
        cutterPlane.SetInputData(InputModel)
        cutterPlane.Update()

        cutterModel = vtk.vtkPolyData()
        cutterModel = cutterPlane.GetOutput()
        surfPlane = vtk.vtkSurfaceReconstructionFilter()
        surfPlane.SetInputData(cutterModel)
        cfPlane = vtk.vtkContourFilter()
        cfPlane.SetInputConnection(surfPlane.GetOutputPort())
        cfPlane.SetValue(0, 0.0)
        reconstructedPlane = vtk.vtkReverseSense()
        reconstructedPlane.SetInputConnection(cfPlane.GetOutputPort())
        reconstructedPlane.ReverseCellsOn()
        reconstructedPlane.ReverseNormalsOn()


        # create a loop defined by the input points
        #PointsPolyData = vtk.vtkPolyData()
        #self.FiducialsToPolyData(fidList, PointsPolyData)

        PointsPolyData = vtk.vtkPolyData()
        #self.errorSim(fidList, PointsPolyData)
        self.FiducialsToPolyData(fidList, PointsPolyData)
        # #****************************************************************************
        # spline model to model back of the chest wall
        spline = vtk.vtkParametricSpline()
        splinePoints = vtk.vtkPoints()
        for i in range(PointsPolyData.GetNumberOfPoints()):
            splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))

        spline.SetPoints(splinePoints)
        spline.ClosedOn()
        parametricSpline = vtk.vtkParametricFunctionSource()
        parametricSpline.SetParametricFunction(spline)
        implictSpline = vtk.vtkImplicitPolyDataDistance()
        implictSpline.SetInput(parametricSpline.GetOutput())

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

        splineTransform = vtk.vtkThinPlateSplineTransform()
        splineTransform.SetSourceLandmarks(projectedPoints)
        splineTransform.SetTargetLandmarks(PointsPolyData.GetPoints())
        splineTransform.SetBasisToR()

        TransformedPlane = vtk.vtkTransformPolyDataFilter()
        TransformedPlane.SetInputConnection(reconstructedPlane.GetOutputPort())
        TransformedPlane.SetTransform(splineTransform)

        planeModel = modelsLogic.AddModel(TransformedPlane.GetOutputPort())
        planeModel.GetDisplayNode().SetVisibility(False)
        planeModel.SetName("transformedPlane")

        loop = vtk.vtkImplicitSelectionLoop()
        loop.SetLoop(PointsPolyData.GetPoints())
        loop.SetNormal(plane.GetNormal())

        implictSplinePlane = vtk.vtkImplicitPolyDataDistance()
        implictSplinePlane.SetInput(TransformedPlane.GetOutput())
    #######################################################################################

        extrude = vtk.vtkLinearExtrusionFilter()
        extrude.SetInputData(TransformedPlane.GetOutput())
        extrude.SetScaleFactor(1)
        extrude.CappingOn()
        normVec = plane.GetNormal()
        if LeftBreast == False:
            extrude.SetVector((normVec[0]), (normVec[1]), (normVec[2]))
        else:
            extrude.SetVector((normVec[0] * -1), (normVec[1] * -1), (normVec[2] * -1))
        extrude.Update()

        implictextrude = vtk.vtkImplicitPolyDataDistance()
        implictextrude.SetInput(extrude.GetOutput())

        # Clip the clipped input model with the loop


        clippedInputWithLoop = vtk.vtkClipPolyData()
        clippedInputWithLoop.SetClipFunction(loop)  # should be loop
        clippedInputWithLoop.SetInputData(InputModel)
        clippedInputWithLoop.SetInsideOut(True)
        clippedInputWithLoop.Update()

        # No use the vtkPolyDataConnectivityFilter to extract the largest region
        connectedInput = vtk.vtkPolyDataConnectivityFilter()
        connectedInput.SetInputConnection(clippedInputWithLoop.GetOutputPort())
        connectedInput.SetExtractionModeToLargestRegion()
        connectedInput.Update()


        # clippedInputWithExtrude = vtk.vtkClipPolyData()
        # clippedInputWithExtrude.SetClipFunction(implictSpline)  # should be loop
        # clippedInputWithExtrude.SetInputData(connectedInput.GetOutput())
        # clippedInputWithExtrude.SetInsideOut(True)
        # clippedInputWithExtrude.Update()


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
        planeModel = modelsLogic.AddModel(clipClosedBreast.GetOutputPort())
        planeModel.GetDisplayNode().SetVisibility(True)
        planeModel.SetName(name)
        planeModel.GetDisplayNode().BackfaceCullingOff()

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

    def cropWithCurve(self, modelNode, fidList, LeftBreast, volume, surfaceArea, modelFids):
        #Check which breast volume is being computed for
        if LeftBreast == True:
            name = "ClosedLeftBreast"
        else:
            name = "ClosedRightBreast"

        modelsLogic = slicer.modules.models.logic()

        PointsPolyData = vtk.vtkPolyData()
        self.errorSim(fidList, PointsPolyData)
        #self.FiducialsToPolyData(fidList, PointsPolyData)

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

        modelsLogic = slicer.modules.models.logic()

        # create a loop defined by the input points
       # PointsPolyData = vtk.vtkPolyData()
        #self.FiducialsToPolyData(fidList, PointsPolyData)

        # modelPoints = vtk.vtkPolyData()
        # self.FiducialsToPolyData(modelFids, modelPoints)

        # #****************************************************************************
        # spline model to model back of the chest wall
        spline = vtk.vtkParametricSpline()
        splinePoints = vtk.vtkPoints()
        # for i in range(modelPoints.GetNumberOfPoints()):
        #     splinePoints.InsertNextPoint(modelPoints.GetPoint(i))
        for i in range(PointsPolyData.GetNumberOfPoints()):
            splinePoints.InsertNextPoint(PointsPolyData.GetPoint(i))

        spline.SetPoints(splinePoints)
        #spline.ClosedOn()
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
        #spline2.ClosedOn()
        functionSource2 = vtk.vtkParametricFunctionSource()
        functionSource2.SetParametricFunction(spline2)

        splineTransform = vtk.vtkThinPlateSplineTransform()
        splineTransform.SetSourceLandmarks(projectedPoints)
        splineTransform.SetTargetLandmarks(PointsPolyData.GetPoints())#(modelPoints.GetPoints())
        splineTransform.SetBasisToR() #Since our points are 3D

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
        clippedInputWithLoop.SetClipFunction(loop)  # should be loop
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
        ########Change this line below
        #************************************************************
        if LeftBreast == True:
            clippedInputhWithPlane.SetInsideOut(False)
        else:
            clippedInputhWithPlane.SetInsideOut(True)
        clippedInputhWithPlane.Update()

        clippedNormals = vtk.vtkPolyDataNormals()
        clippedNormals.SetInputConnection(clipped.GetOutputPort())
        clippedNormals.ComputePointNormalsOn()
        clippedNormals.ConsistencyOn()
        #If line 509 is changed to true then line 517 must also be commented out
        ####change this line below
        #*******************************************
        #clippedNormals.FlipNormalsOn()
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

        return finalModel

    ############################################################
    #All should be removed only to get image with cropped out breast and curved surface
        # CenterOfMass = vtk.vtkCenterOfMass()
        # CenterOfMass.SetInputData(PointsPolyData)
        # CenterOfMass.SetUseScalarsAsWeights(False)
        # CenterOfMass.Update()
        # center = CenterOfMass.GetCenter()
        #
        #
        # sphere1 = vtk.vtkSphere()
        # sphere1.SetCenter(center)
        # sphere1.SetRadius(87)
        #
        # sphere2 = vtk.vtkSphere()
        # sphere2.SetCenter(center)
        # sphere2.SetRadius(83)
        #
        # clippedInputWithLoop2 = vtk.vtkClipPolyData()
        # clippedInputWithLoop2.SetClipFunction(sphere2)  # should be loop
        # clippedInputWithLoop2.SetInputData(InputModel)
        # clippedInputWithLoop2.SetInsideOut(False)
        # clippedInputWithLoop2.Update()
        #
        # clippedInputWithLoop3 = vtk.vtkClipPolyData()
        # clippedInputWithLoop3.SetClipFunction(sphere1)  # should be loop
        # clippedInputWithLoop3.SetInputData(TransformedPlane.GetOutput())
        ######Change here as well################################
        # clippedInputWithLoop3.SetInsideOut(True)
        # clippedInputWithLoop3.Update()
        #
        # appendClosedBreast2 = vtk.vtkAppendPolyData()
        # appendClosedBreast2.AddInputData(clippedInputWithLoop2.GetOutput())  # need to change these variables
        # appendClosedBreast2.AddInputData(clippedInputWithLoop3.GetOutput())
        # appendClosedBreast2.Update()
        #
        # finalModel = modelsLogic.AddModel(clippedInputWithLoop2.GetOutput())
        # finalModel.GetDisplayNode().SetVisibility(True)
        # finalModel.SetName("surface")
        # finalModel.GetDisplayNode().BackfaceCullingOff()
        #
        # finalModel = modelsLogic.AddModel(clippedInputWithLoop3.GetOutput())
        # finalModel.GetDisplayNode().SetVisibility(True)
        # finalModel.SetName("surface2")
        # finalModel.GetDisplayNode().BackfaceCullingOff()

    def errorSim(self, fidList, polyData):
        PointsPolyData = vtk.vtkPolyData()
        self.FiducialsToPolyData(fidList, PointsPolyData)
        NumberOfPoints = PointsPolyData.GetNumberOfPoints()
        points = vtk.vtkPoints()
        for i in range(NumberOfPoints):
            p = PointsPolyData.GetPoint(i)
            ran1 = uniform(0,10) -5
            ran2 = uniform(0,10) -5
            ran3 = uniform(0,10) -5
            p1 = p[0] + ran1
            p2 = p[1] + ran2
            p3 = p[2] + ran3
            points.InsertNextPoint([p1,p2,p3])
            tempPolyData = vtk.vtkPolyData()
            tempPolyData.SetPoints(points)

            vertex = vtk.vtkVertexGlyphFilter()
            vertex.SetInputData(tempPolyData)
            vertex.Update()

            polyData.ShallowCopy(vertex.GetOutput())


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
        import vtkSegmentationCorePython as vtkSegmentationCore

        # segmentationNode = slicer.mrmlScene.GetNthNodeByClass(1, 'vtkMRMLSegmentationNode')
        # if segmentationNode == None:
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
        return volume




    def run(self, inputModel, fidList, LeftBreast, volume, surfaceArea,  planeType, modelFids):
        """
        Run the actual algorithm
        """
        #If flat cut plane is chosen
        if planeType == True:
            self.cropWithPlane(inputModel, fidList, LeftBreast, volume, surfaceArea)
            self.AddVolumeNode()
            return 0
        #else curved cut plane is chosen
        else:
            modelNode = self.cropWithCurve(inputModel, fidList, LeftBreast, volume, surfaceArea, modelFids)
            self.AddVolumeNode()
            volumeVal = self.AddSegmentNode(modelNode, volume)
            return volumeVal

        logging.info('Processing completed')



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