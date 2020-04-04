# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:15:44 2020

@author: Arsenic
"""

from paraview.simple import *
import numpy as np
from math import pi, cos, sin

band = 1484
grain = 1217
origin = np.genfromtxt("K:\\Marie\\718\\origin.txt")  # read origin from file
normal = [0.75, 0.62, -0.23]  # read normal from ext file

def deg_to_rad(deg):
    rad = deg*pi/180
    return rad

def rad_to_deg(rad):
    deg = rad*180/pi
    return deg

def rot_matrix(eulers):
    "macro base in crystal frame"
    "RD = rotated_matrix[0,:]"
    "TD = rotated_matrix[1,:]"
    "ND = rotated_matrix[2,:]"
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
    T = np.array([[0,-1,0],[0,0,1],[1,0,0]])
    rotation_matrix3D = np.dot(rotation_matrix2D,T)
    return rotation_matrix3D

def create_threshold(grainID):
    # find source
    a718_final6xdmf = FindSource('718_final6.xdmf')
    # create a new 'Threshold'
    threshold5 = Threshold(Input=a718_final6xdmf)
    threshold5.Scalars = ['CELLS', 'FeatureIdsClean']
    # Properties of threshold5
    threshold5.ThresholdRange = [grainID, grainID]
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1387, 548]
    # show data in view
    threshold5Display = Show(threshold5, renderView1)
    # get color transfer function/color map for 'FeatureIdsClean'
    featureIdsCleanLUT = GetColorTransferFunction('FeatureIdsClean')
    # get opacity transfer function/opacity map for 'FeatureIdsClean'
    featureIdsCleanPWF = GetOpacityTransferFunction('FeatureIdsClean')

# trace defaults for the display properties.
    threshold5Display.Representation = 'Surface'
    threshold5Display.ColorArrayName = ['CELLS', 'FeatureIdsClean']
    threshold5Display.LookupTable = featureIdsCleanLUT
    threshold5Display.EdgeColor = [1.0, 1.0, 1.0]
    threshold5Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold5Display.SelectOrientationVectors = 'IPFColorX'
    threshold5Display.ScaleFactor = 21.700000000000003
    threshold5Display.SelectScaleArray = 'FeatureIdsClean'
    threshold5Display.GlyphType = 'Arrow'
    threshold5Display.GlyphTableIndexArray = 'FeatureIdsClean'
    threshold5Display.GaussianRadius = 1.085
    threshold5Display.SetScaleArray = [None, '']
    threshold5Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold5Display.OpacityArray = [None, '']
    threshold5Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold5Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold5Display.SelectionCellLabelFontFile = ''
    threshold5Display.SelectionPointLabelFontFile = ''
    threshold5Display.PolarAxes = 'PolarAxesRepresentation'
    threshold5Display.ScalarOpacityFunction = featureIdsCleanPWF
    threshold5Display.ScalarOpacityUnitDistance = 2.9060731563875963

    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    threshold5Display.DataAxesGrid.XTitleFontFile = ''
    threshold5Display.DataAxesGrid.YTitleFontFile = ''
    threshold5Display.DataAxesGrid.ZTitleFontFile = ''
    threshold5Display.DataAxesGrid.XLabelFontFile = ''
    threshold5Display.DataAxesGrid.YLabelFontFile = ''
    threshold5Display.DataAxesGrid.ZLabelFontFile = ''

    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    threshold5Display.PolarAxes.PolarAxisTitleFontFile = ''
    threshold5Display.PolarAxes.PolarAxisLabelFontFile = ''
    threshold5Display.PolarAxes.LastRadialAxisTextFontFile = ''
    threshold5Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

    # hide data in view
    Hide(a718_final6xdmf, renderView1)
    # show color bar/color legend
    threshold5Display.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # set scalar coloring
    ColorBy(threshold5Display, ('CELLS', 'IPFColorX', 'Magnitude'))
    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(featureIdsCleanLUT, renderView1)
    # rescale color and/or opacity maps used to include current data range
    threshold5Display.RescaleTransferFunctionToDataRange(True, False)
    # show color bar/color legend
    threshold5Display.SetScalarBarVisibility(renderView1, True)
        # get color transfer function/color map for 'IPFColorX'
    iPFColorXLUT = GetColorTransferFunction('IPFColorX')
    # get opacity transfer function/opacity map for 'IPFColorX'
    iPFColorXPWF = GetOpacityTransferFunction('IPFColorX')
    # Properties modified on threshold5Display
    threshold5Display.MapScalars = 0
    # Properties modified on threshold5Display
    threshold5Display.Opacity = 0.3
    # rename source object
    RenameSource('grain %s' %grainID, threshold5)
    #### saving camera placements for all active views
    # current camera placement for renderView1
    renderView1.CameraPosition = [266.2080813646024, -259.50780354197593, 574.7555833106233]
    renderView1.CameraFocalPoint = [172.50000000000006, 63.50000000000004, 137.50000000000006]
    renderView1.CameraViewUp = [0.07608517732440345, -0.7941416604379515, -0.6029511331346574]
    renderView1.CameraParallelScale = 142.77517291181965

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
    
def draw_slice(origin, normal, bandID, grainID):
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

    # get color transfer function/color map for 'FeatureIdsClean'
    featureIdsCleanLUT = GetColorTransferFunction('FeatureIdsClean')

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = ['CELLS', 'FeatureIdsClean']
    slice1Display.LookupTable = featureIdsCleanLUT
    slice1Display.EdgeColor = [1.0, 1.0, 1.0]
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'IPFColorX'
    slice1Display.ScaleFactor = 20.1
    slice1Display.SelectScaleArray = 'FeatureIdsClean'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'FeatureIdsClean'
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
    
    # get opacity transfer function/opacity map for 'FeatureIdsClean'
    featureIdsCleanPWF = GetOpacityTransferFunction('FeatureIdsClean')
    
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
    

create_threshold(grain)
draw_slice(origin,normal,band,grain)
