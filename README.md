# BreastReconstruction
*Tool for breast reconstruction surgery planning*

A 3D Slicer module to measure the breast volume difference given a 3D surface scan and breast boundaries as input.
The module can be used to measure the breast volume difference between the right and left breast from the same surface scan or between the same breast from multiple scans. 

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
  
## Workflow for breast volume difference between left and right breast 
1. Open the BreastRreconstruction Module from the modules drop-down menu 
2. Import the 3D surface scan for which you would like to measure the breast volume difference (drag scan into 3D Slicer) 
3. Either:
    1. Import predefined left and right breast fiduicals for the patient 
    2. Create new breast boundary fiduicals see [Creating Breast Boundary](#breastBond)
3. Select the surface scan as the input model
4. Select the breast boundary fiducials as the left and right fiducials
4. Click the *Apply* Button
*The module should now display the left breast volume, right breast volume and the breast volume difference. The module will also create to new models one called ClosedRightBreast and the other ClosedleftBreast*
5. Verify the module has performed correctly
    1. Open the Models module
    2. Inspect the ClosedRightBreast and ClosedLeftBreast, they should repersent the breasts clipped from the model with a curved posterior wall. Below is an image of the module performing correctly with the sample mannequin scan. [picture alt]()

## 3D Surface Scans 

## Creating Breast Boundary <a name="breastBond"></a>


