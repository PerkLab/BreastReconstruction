cmake_minimum_required(VERSION 2.8.9)

project(BreastReconstruction)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/BreastReconstruction")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Rachael House (Perk Lab)")
set(EXTENSION_DESCRIPTION "Extension to calculate the volume of the breast")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/BreastReconstruction.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/BreastReconstruction/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BreastReconstruction)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
