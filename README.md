# BreastReconstruction
*Tool for breast reconstruction surgery planning*

A 3D Slicer module to measure the breast volume difference given a 3D surface scan and breast boundaries markups as input.

## Installation
1. Clone or download the Breast Reconstruction Module [here](https://github.com/PerkLab/BreastReconstruction)
2. Install 3D slicer from the [download page](https://download.slicer.org/)
3. Install the Breast Reconstruction Module: 
    - Open 3D Slicer
    - Click *Edit* in the top left corner 
    - Select *Application Settings* from the drop-down list
    - Select *Modules* from list on the left-hand side of the panel
    - Click *Add* under Paths on the right-hand side of the panel
    - Select path to where the Breast Reconstruction module is downloaded/cloned on your local machine (select path to depth of the .py files)
    - Select *OK* bottom of the panel
    - Restart 3D Silcer 
    
*The module should now be installed and visible under the Modules drop-down menu under the category BreastSurgery.* 
  
## Workflow
1. Open the *BreastRreconstruction* Module from the modules drop-down menu. Below is an screen shot of the modules user interface.

![](https://github.com/PerkLab/BreastReconstruction/blob/master/data/ExampleScreenshots/userInterface.PNG "User Interface")

2. Import the 3D surface scan for which you would like to measure the breast volume difference (drag scan into 3D Slicer) 
3. Either:
    1. If this is the intital scan, create the new breast boundary fiduicals:
        1. Select *Create new MarkupFiducials* from the dropdown menu assosiated with the left or right fiduicals 
        2. Place fiduicals around the circumference of the breast starting from middle of the chest
        3. Repeat for left and right breast
    2. Register to the intital scan and import fiducials 
        1. Import intital scan and use it as input to optional 
        2. Input left and right fiducials and put them in corresponding input
        
        *When creating the breast boundary it is important that is remains constant for each patient, the boundary should be marked once and the fiduicals saved for other measurements. The breast boundary should resemble the breast footprint, this is done by having the breast reconstruction surgeon mark the footprint on the patient before the initial scan.*
4. Select the surface scan as the input model
5. Click the *Apply* Button

*The module should now display the left breast volume, right breast volume and the breast volume difference. The module will also create two new models one called ClosedRightBreast and the other ClosedleftBreast*

6. Verify the module has performed correctly
    1. Open the *Models* module
    2. Inspect the ClosedRightBreast and ClosedLeftBreast, they should represent the breasts clipped from the model with a curved posterior wall. Below are an images of the module performing correctly with the sample mannequin scan. 
    
    ![](https://github.com/PerkLab/BreastReconstruction/blob/master/data/ExampleScreenshots/manequinBreastsFront.PNG "Front View")![](https://github.com/PerkLab/BreastReconstruction/blob/master/data/ExampleScreenshots/manequinBreastsSide.PNG "Side View")
    
    3. If the module does not produce the correct results please follow the instructions [here](#error)

## Error Handling <a name="error"></a>
If after inspection the closed right and/or left breast models are produced please consult below.

There are two options that can be used together or seperate for both breasts to correct the meshes.
1. Reverse the normal of the breast. An example of a breast where the normal was reversed to obtain the correct results is below.  
![](https://github.com/PerkLab/BreastReconstruction/blob/master/data/ExampleScreenshots/reverseNormal.PNG "Reverse Normal Example")
2. Set the clip function to inside out. An example of a breast where the the clip function was set to inside out to obtain the correct results is below.  
![](https://github.com/PerkLab/BreastReconstruction/blob/master/data/ExampleScreenshots/insideOutClip.PNG "Inside Out Clip Function Example")

## 3D Surface Scans 
The 3D surface scans for this project were captured using the [Artec Eva](https://www.artec3d.com/). The surface scans are processed using the Artec Studio software. The scans are saved as .obj files and the scans textures are saved as .jpeg files. 

## Example Scan
A mannequin was scanned using the Artec Eva, the scan is stored [here](https://github.com/PerkLab/BreastReconstruction/tree/master/data/Example3Dscans). Also stored at this location is an example of a left and right breast fiduical set for the mannequin.



