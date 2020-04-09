# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:15:44 2020

@author: Arsenic
"""

from paraview.simple import *
import numpy as np
from math import pi, cos, sin

dataset_name = '718_final6.xdmf'
grain_names = 'FeatureIdsClean'
display = 'IPFColorX'
band_list = [1913,448]
step = 0
plotWholeTRD = True  # make true if you want to plot the whole TRD containing the band for more realistic shape


# microstructure information
grainHosts = np.load("K:\\Marie\\718\\grainHost.npy")
parentTRD = np.genfromtxt("K:\\Marie\\718\\ParentIds.txt")  # load parent IDs 
feat2Dto3D = np.genfromtxt("K:\\Marie\\718\\feat2Dto3D.txt")  # correspondance between 2D and 3D features
bandendpoints = np.load("K:\\Marie\\718\\bandEndpoints.npy")
hklPlanes = np.load("K:\\Marie\\718\\hklPlane.npy")
eulers = np.genfromtxt("K:\\Marie\\718\\grainEulers.txt")

def deg_to_rad(deg):
    rad = deg*pi/180
    return rad

def rad_to_deg(rad):
    deg = rad*180/pi
    return deg

def rot_matrix(eulers):
    "ANGLES HAVE TO BE IN RAD!"
    "macro base in crystal frame"
    "RD = rotated_matrix[0,:]"
    "TD = rotated_matrix[1,:]"
    "ND = rotated_matrix[2,:]"
    if eulers.max() > 6.30:
        eulers = deg_to_rad(eulers)
    rotated_matrix = np.zeros((3,3))
    rotated_matrix[0,0] = cos(eulers[0])*cos(eulers[2])-sin(eulers[0])*sin(eulers[2])*cos(eulers[1])
    rotated_matrix[0,1] = -cos(eulers[0])*sin(eulers[2])-sin(eulers[0])*cos(eulers[2])*cos(eulers[1])
    rotated_matrix[0,2] = sin(eulers[0])*sin(eulers[1])
    rotated_matrix[1,0] = sin(eulers[0])*cos(eulers[2])+cos(eulers[0])*sin(eulers[2])*cos(eulers[1])
    rotated_matrix[1,1] = -sin(eulers[0])*sin(eulers[2])+cos(eulers[0])*cos(eulers[2])*cos(eulers[1])
    rotated_matrix[1,2] = -cos(eulers[0])*sin(eulers[1])
    rotated_matrix[2,0] = sin(eulers[2])*sin(eulers[1])
    rotated_matrix[2,1] = cos(eulers[2])*sin(eulers[1])
    rotated_matrix[2,2] = cos(eulers[1])
    return rotated_matrix

def gen_plane(plane_number):
    if plane_number == 1:
        plane = np.array([1,1,1])
    if plane_number == 2:
        plane = np.array([-1,1,1])
    if plane_number == 3:
        plane = np.array([1,-1,1])
    if plane_number == 4:
        plane = np.array([1,1,-1])
    if plane_number == 0:
        plane = np.array([0,0,0])
    return plane

def get_plane_normal(plane_hkl,rotated_matrix):
    # plane_normal = np.zeros((3),dtype=float)
    plane_normal = np.dot(rotated_matrix,plane_hkl)
    # plane_normal[0] = (np.cross(plane_hkl,rotated_matrix[2,:]),rotated_matrix[0,:])
    # plane_normal[1] = np.dot(np.cross(plane_hkl,rotated_matrix[2,:]),rotated_matrix[1,:])
    # plane_normal[2] = np.dot(np.cross(plane_hkl,rotated_matrix[2,:]),rotated_matrix[2,:])
    return plane_normal

def rotateVector(p, q):
	## returns vector rotated by a quaternion
	## p: the vector to rotate
	## q: the quaternion
	## quaternion of the form (x, y, z, w), this is the DREAM.3D ordering
	## uses p' = qpq*
	## details of math at https://www.essentialmath.com/GDC2013/GDC13_quaternions_final.pdf
	v = q[:3]
	w = q[3]
	pRotated = p + (2 * w) * (np.cross(v, p)) + (2 * (np.cross(v, (np.cross(v,p)))))
	return pRotated

def find_normal_3D(normal_2D):
    normal_3D = np.zeros_like(normal_2D)
    normal_3D[0] = -normal_2D[1]
    normal_3D[1] = normal_2D[2]
    normal_3D[2] = normal_2D[0]
    return normal_3D

def rotate_to_3D(rotation_matrix2D):
    # T = np.array([[0,0,1],[-1,0,0],[0,-1,0]])  #
    T = np.array([[0,-1,0],[0,0,-1],[1,0,0]])
    rotation_matrix3D = np.dot(T,rotation_matrix2D)
    return rotation_matrix3D

def norm(v):
    if np.sqrt(np.sum(v**2)) == 0:
        normalized_v = v
    else:
        normalized_v = v / np.sqrt(np.sum(v**2))
    return normalized_v

def create_threshold(grainID):
    # find source
    xdmf_file = FindSource(dataset_name)
    # create a new 'Threshold'
    threshold = Threshold(Input=xdmf_file)
    threshold.Scalars = ['CELLS', grain_names]
    # Properties of threshold
    threshold.ThresholdRange = [grainID, grainID]
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1387, 548]
    # show data in view
    thresholdDisplay = Show(threshold, renderView1)
    # get color transfer function/color map for grain_names
    featureIdsCleanLUT = GetColorTransferFunction(grain_names)
    # get opacity transfer function/opacity map for grain_names
    featureIdsCleanPWF = GetOpacityTransferFunction(grain_names)

# trace defaults for the display properties.
    thresholdDisplay.Representation = 'Surface'
    thresholdDisplay.ColorArrayName = ['CELLS', grain_names]
    thresholdDisplay.LookupTable = featureIdsCleanLUT
    thresholdDisplay.EdgeColor = [1.0, 1.0, 1.0]
    thresholdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    thresholdDisplay.SelectOrientationVectors = display
    thresholdDisplay.ScaleFactor = 21.700000000000003
    thresholdDisplay.SelectScaleArray = grain_names
    thresholdDisplay.GlyphType = 'Arrow'
    thresholdDisplay.GlyphTableIndexArray = grain_names
    thresholdDisplay.GaussianRadius = 1.085
    thresholdDisplay.SetScaleArray = [None, '']
    thresholdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    thresholdDisplay.OpacityArray = [None, '']
    thresholdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    thresholdDisplay.DataAxesGrid = 'GridAxesRepresentation'
    thresholdDisplay.SelectionCellLabelFontFile = ''
    thresholdDisplay.SelectionPointLabelFontFile = ''
    thresholdDisplay.PolarAxes = 'PolarAxesRepresentation'
    thresholdDisplay.ScalarOpacityFunction = featureIdsCleanPWF
    thresholdDisplay.ScalarOpacityUnitDistance = 2.9060731563875963

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    thresholdDisplay.DataAxesGrid.XTitleFontFile = ''
    thresholdDisplay.DataAxesGrid.YTitleFontFile = ''
    thresholdDisplay.DataAxesGrid.ZTitleFontFile = ''
    thresholdDisplay.DataAxesGrid.XLabelFontFile = ''
    thresholdDisplay.DataAxesGrid.YLabelFontFile = ''
    thresholdDisplay.DataAxesGrid.ZLabelFontFile = ''

    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    thresholdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    thresholdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    thresholdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    thresholdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

    # hide data in view
    Hide((xdmf_file), renderView1)
    # show color bar/color legend
    thresholdDisplay.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # set scalar coloring
    ColorBy(thresholdDisplay, ('CELLS', display, 'Magnitude'))
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(featureIdsCleanLUT, renderView1)
    # rescale color and/or opacity maps used to include current data range
    thresholdDisplay.RescaleTransferFunctionToDataRange(True, False)
    # show color bar/color legend
    thresholdDisplay.SetScalarBarVisibility(renderView1, True)
        # get color transfer function/color map for display
    displayLUT = GetColorTransferFunction(display)
    # get opacity transfer function/opacity map for display
    displayPWF = GetOpacityTransferFunction(display)
    # Properties modified on thresholdDisplay
    thresholdDisplay.MapScalars = 0
    # Properties modified on thresholdDisplay
    thresholdDisplay.Opacity = 0.3
    # rename source object
    RenameSource('grain %s' %grainID, threshold)
    #### saving camera placements for all active views
    # current camera placement for renderView1
    renderView1.CameraPosition = [266.2080813646024, -259.50780354197593, 574.7555833106233]
    renderView1.CameraFocalPoint = [172.50000000000006, 63.50000000000004, 137.50000000000006]
    renderView1.CameraViewUp = [0.07608517732440345, -0.7941416604379515, -0.6029511331346574]
    renderView1.CameraParallelScale = 142.77517291181965

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
    
def draw_slice(origin, normal, bandID,grainID):
    paraview.simple._DisableFirstRenderCameraReset()
    grain = FindSource('grain %s' %grainID)
    # create a new 'Slice'
    slice1 = Slice(Input=grain)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # set active source
    SetActiveSource(slice1)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [790, 548]

    # show data in view
    slice1Display = Show(slice1, renderView1)

    # get color transfer function/color map for grain_names
    featureIdsCleanLUT = GetColorTransferFunction(grain_names)

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = ['CELLS', grain_names]
    slice1Display.LookupTable = featureIdsCleanLUT
    slice1Display.EdgeColor = [1.0, 1.0, 1.0]
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = display
    slice1Display.ScaleFactor = 20.1
    slice1Display.SelectScaleArray = grain_names
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = grain_names
    slice1Display.GaussianRadius = 1.0050000000000001
    slice1Display.SetScaleArray = [None, '']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = [None, '']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.SelectionCellLabelFontFile = ''
    slice1Display.SelectionPointLabelFontFile = ''
    slice1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    slice1Display.DataAxesGrid.XTitleFontFile = ''
    slice1Display.DataAxesGrid.YTitleFontFile = ''
    slice1Display.DataAxesGrid.ZTitleFontFile = ''
    slice1Display.DataAxesGrid.XLabelFontFile = ''
    slice1Display.DataAxesGrid.YLabelFontFile = ''
    slice1Display.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
    slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
    slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
    slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)
    
    # get opacity transfer function/opacity map for grain_names
    featureIdsCleanPWF = GetOpacityTransferFunction(grain_names)
    
    # turn off scalar coloring
    ColorBy(slice1Display, None)
    
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(featureIdsCleanLUT, renderView1)
    
    # change solid color
    slice1Display.DiffuseColor = [0.0, 0.3333333333333333, 0.4980392156862745]

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = origin  # 
    slice1.SliceType.Normal = normal  #

    # show data in view
    slice1Display = Show(slice1, renderView1)
    # rename source object
    RenameSource('band %s' %bandID, slice1)

    # update the view to ensure updated data information
    renderView1.Update()
    #### saving camera placements for all active views
    # current camera placement for renderView1
    renderView1.CameraPosition = [436.9555297314597, -177.82337347040382, 557.1832055368379]
    renderView1.CameraFocalPoint = [172.5, 63.5, 137.5]
    renderView1.CameraViewUp = [-0.104225946984368, -0.889146770559922, -0.44559507669864734]
    renderView1.CameraParallelScale = 142.77517291181965


def plot_all_bands(band_list):
    for band in band_list:
        grainHost2D = int(grainHosts[step,band])  # grain of band in 2D
        grainHost3D = int(feat2Dto3D[grainHost2D])  # corresponding 3D grain
        grain_list = []
        for grain in range(parentTRD.shape[0]):
            if parentTRD[grain] == parentTRD[grainHost3D]:
                grain_list.append(grain)
        origin_2D = bandendpoints[step,band,:2]  # pick any of the 2 points
        origin_3D = [origin_2D[0],0,origin_2D[1]]  # y=0 is the DIC surface
        plane_hkl = gen_plane(hklPlanes[step,band])  # get plane ID
        eulers2D = eulers[grainHost2D]
        rotMatrix2D = rot_matrix(eulers2D)
        normal_2D = np.dot(rotMatrix2D,plane_hkl)
        normal_3D = norm(np.array([-normal_2D[1],normal_2D[2],normal_2D[0]]))
        if plotWholeTRD == True:
            for grain in grain_list:
                create_threshold(grain)
        if plotWholeTRD == False:
            create_threshold(grainHost3D)  # create threshold corresponding to the grain of interest
        draw_slice(origin_3D,normal_3D,band,grainHost3D)
# draw_slice(origin_3D,normal_3D,bandID_,grainID_)
plot_all_bands(band_list)