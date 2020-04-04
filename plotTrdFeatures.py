# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 18:02:05 2020

@author: Pollock
"""


from paraview.simple import *
import numpy as np

grainID = 1217

# load parent feature IDs
parentTRD = np.genfromtxt("K:\\Marie\\718\\ParentIds.txt")
grains = [1217]
for grain in range(parentTRD.shape[0]):
    if parentTRD[grain] == parentTRD[grainID]:
        grains.append(grain)

def plot_feature_trd(grainID):
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

for grain in grains:
    plot_feature_trd(grain)