# state file generated using paraview version 5.12.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

if len(sys.argv) == 2:
        dim = int(sys.argv[1])
else:
        dim = 128

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [432, 778]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'Grid Axes 3D Actor'
renderView2.OrientationAxesVisibility = 0
renderView2.CenterOfRotation = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [2035.5000274181366, 2035.5000274181366, 26195.030829693267]
renderView2.CameraFocalPoint = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 5512.570271480759
renderView2.LegendGrid = 'Legend Grid Actor'
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.ViewSize = [434, 778]
renderView4.AxesGrid = 'Grid Axes 3D Actor'
renderView4.OrientationAxesVisibility = 0
renderView4.CenterOfRotation = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView4.StereoType = 'Crystal Eyes'
renderView4.CameraPosition = [1783.8760241992409, 705.2113322092382, 26160.020212472777]
renderView4.CameraFocalPoint = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView4.CameraViewUp = [-0.0001626974090163947, 0.9987096067050545, 0.050784791075561234]
renderView4.CameraViewAngle = 23.70600414078675
renderView4.CameraFocalDisk = 1.0
renderView4.CameraParallelScale = 6976.169714738533
renderView4.LegendGrid = 'Legend Grid Actor'
renderView4.BackEnd = 'OSPRay raycaster'
renderView4.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView5 = CreateView('RenderView')
renderView5.ViewSize = [432, 778]
renderView5.AxesGrid = 'Grid Axes 3D Actor'
renderView5.OrientationAxesVisibility = 0
renderView5.CenterOfRotation = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView5.StereoType = 'Crystal Eyes'
renderView5.CameraPosition = [2035.5000274181366, 2035.5000274181366, 26195.030829693267]
renderView5.CameraFocalPoint = [2035.5000274181366, 2035.5000274181366, 0.0]
renderView5.CameraViewAngle = 23.39544513457557
renderView5.CameraFocalDisk = 1.0
renderView5.CameraParallelScale = 6976.169714738533
renderView5.LegendGrid = 'Legend Grid Actor'
renderView5.BackEnd = 'OSPRay raycaster'
renderView5.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.333333)
layout1.AssignView(1, renderView4)
layout1.SplitHorizontal(2, 0.500000)
layout1.AssignView(5, renderView2)
layout1.AssignView(6, renderView5)
layout1.SetSize(1300, 778)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView4)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
backpackpvd = PVDReader(registrationName='backpack.pvd', FileName='backpack.pvd')
backpackpvd.PointArrays = ['ImageFile']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=backpackpvd)
resampleToImage1.SamplingDimensions = [dim, dim, dim]
resampleToImage1.SamplingBounds = [0.0, 511.0, 0.0, 511.0, 0.0, 372.0]

ghostCellsGenerator1 = GhostCellsGenerator(registrationName='GhostCellsGenerator1', Input=resampleToImage1)
ghostCellsGenerator1.MinimumNumberOfGhosts = 1

tTKArrayPreconditioning = TTKArrayPreconditioning(registrationName='TTKArrayPreconditioning5', Input=ghostCellsGenerator1)
tTKArrayPreconditioning.PointDataArrays = ['ImageFile']
tTKArrayPreconditioning.GlobalOrderArray = 1

# create a new 'TTK PersistenceDiagram'
tTKPersistenceDiagram1 = TTKPersistenceDiagram(registrationName='TTKPersistenceDiagram1', Input=tTKArrayPreconditioning)
tTKPersistenceDiagram1.ScalarField = ['POINTS', 'ImageFile']
tTKPersistenceDiagram1.InputOffsetField = ['POINTS', 'ImageFile']

SaveData('diagram.pvtu', proxy=tTKPersistenceDiagram1)

# create a new 'Threshold'
threshold5 = Threshold(registrationName='Threshold5', Input=tTKPersistenceDiagram1)
threshold5.Scalars = ['CELLS', 'PairType']
threshold5.LowerThreshold = 1.0
threshold5.UpperThreshold = 1.0

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=tTKPersistenceDiagram1)
threshold2.Scalars = ['CELLS', 'PairType']
threshold2.LowerThreshold = 0.00025549999554641545
threshold2.UpperThreshold = 14881302.0

# create a new 'Threshold'
threshold3 = Threshold(registrationName='Threshold3', Input=threshold2)
threshold3.Scalars = ['CELLS', 'PairType']
threshold3.LowerThreshold = 2.0
threshold3.UpperThreshold = 2.0

# create a new 'Threshold'
threshold4 = Threshold(registrationName='Threshold4', Input=tTKPersistenceDiagram1)
threshold4.Scalars = ['CELLS', 'PairType']
threshold4.LowerThreshold = 0
threshold4.UpperThreshold = 0

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=tTKPersistenceDiagram1)
threshold1.Scalars = ['CELLS', 'PairIdentifier']
threshold1.LowerThreshold = -1.0
threshold1.UpperThreshold = -1.0

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=threshold1)

# create a new 'Tube'
tube1 = Tube(registrationName='Tube1', Input=extractSurface1)
tube1.Scalars = ['POINTS', '']
tube1.Vectors = ['POINTS', 'Coordinates']
tube1.Radius = 10.0

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from tube1
tube1Display = Show(tube1, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = ['POINTS', '']
tube1Display.DiffuseColor = [0.0, 0.0, 0.0]
tube1Display.SelectTCoordArray = 'None'
tube1Display.SelectNormalArray = 'TubeNormals'
tube1Display.SelectTangentArray = 'None'
tube1Display.OSPRayScaleArray = 'Coordinates'
tube1Display.OSPRayScaleFunction = 'Piecewise Function'
tube1Display.Assembly = ''
tube1Display.SelectOrientationVectors = 'Coordinates'
tube1Display.ScaleFactor = 0.1371263176202774
tube1Display.SelectScaleArray = 'Coordinates'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'Coordinates'
tube1Display.GaussianRadius = 0.00685631588101387
tube1Display.SetScaleArray = ['POINTS', 'Coordinates']
tube1Display.ScaleTransferFunction = 'Piecewise Function'
tube1Display.OpacityArray = ['POINTS', 'Coordinates']
tube1Display.OpacityTransferFunction = 'Piecewise Function'
tube1Display.DataAxesGrid = 'Grid Axes Representation'
tube1Display.PolarAxes = 'Polar Axes Representation'
tube1Display.SelectInputVectors = ['POINTS', 'Coordinates']
tube1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
tube1Display.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
tube1Display.ScaleTransferFunction.Points = [66.0, 0.0, 0.5, 0.0, 1981.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
tube1Display.OpacityTransferFunction.Points = [66.0, 0.0, 0.5, 0.0, 1981.0, 1.0, 0.5, 0.0]

# show data from threshold3
threshold3Display = Show(threshold3, renderView2, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'PairType'
pairTypeTF2D = GetTransferFunction2D('PairType')

# get color transfer function/color map for 'PairType'
pairTypeLUT = GetColorTransferFunction('PairType')
pairTypeLUT.TransferFunction2D = pairTypeTF2D
pairTypeLUT.RGBPoints = [-1.0, 0.929412, 1.0, 1.0, -0.5499267578125, 0.439216, 0.611765, 0.729412, -0.099853515625, 0.235294, 0.333333, 0.501961, 0.3502197265624998, 0.066667, 0.078431, 0.2, 0.5027740250363064, 0.2, 0.1, 0.2, 0.6444497194810088, 0.4, 0.039216, 0.058824, 1.025349093017578, 0.890196, 0.411765, 0.019608, 1.6254347473144533, 0.968627, 0.905882, 0.486275, 2.00048828125, 1.0, 1.0, 0.7]
pairTypeLUT.ColorSpace = 'Lab'
pairTypeLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'PairType'
pairTypePWF = GetOpacityTransferFunction('PairType')
pairTypePWF.Points = [-1.0, 0.0, 0.5, 0.0, 2.00048828125, 1.0, 0.5, 0.0]
pairTypePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
threshold3Display.Representation = 'Surface'
threshold3Display.ColorArrayName = ['CELLS', 'PairType']
threshold3Display.LookupTable = pairTypeLUT
threshold3Display.SelectTCoordArray = 'None'
threshold3Display.SelectNormalArray = 'None'
threshold3Display.SelectTangentArray = 'None'
threshold3Display.OSPRayScaleArray = 'Coordinates'
threshold3Display.OSPRayScaleFunction = 'Piecewise Function'
threshold3Display.Assembly = ''
threshold3Display.SelectOrientationVectors = 'Coordinates'
threshold3Display.ScaleFactor = 0.10448978841304779
threshold3Display.SelectScaleArray = 'Coordinates'
threshold3Display.GlyphType = 'Arrow'
threshold3Display.GlyphTableIndexArray = 'Coordinates'
threshold3Display.GaussianRadius = 0.005224489420652389
threshold3Display.SetScaleArray = ['POINTS', 'Coordinates']
threshold3Display.ScaleTransferFunction = 'Piecewise Function'
threshold3Display.OpacityArray = ['POINTS', 'Coordinates']
threshold3Display.OpacityTransferFunction = 'Piecewise Function'
threshold3Display.DataAxesGrid = 'Grid Axes Representation'
threshold3Display.PolarAxes = 'Polar Axes Representation'
threshold3Display.ScalarOpacityFunction = pairTypePWF
threshold3Display.ScalarOpacityUnitDistance = 0.010490967301258645
threshold3Display.OpacityArrayName = ['POINTS', 'Coordinates']
threshold3Display.SelectInputVectors = ['POINTS', 'Coordinates']
threshold3Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold3Display.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2047.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2047.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from tube1
tube1Display_1 = Show(tube1, renderView4, 'GeometryRepresentation')

# trace defaults for the display properties.
tube1Display_1.Representation = 'Surface'
tube1Display_1.ColorArrayName = ['POINTS', '']
tube1Display_1.DiffuseColor = [0.0, 0.0, 0.0]
tube1Display_1.PointSize = 1.0
tube1Display_1.SelectTCoordArray = 'None'
tube1Display_1.SelectNormalArray = 'TubeNormals'
tube1Display_1.SelectTangentArray = 'None'
tube1Display_1.OSPRayScaleArray = 'Coordinates'
tube1Display_1.OSPRayScaleFunction = 'Piecewise Function'
tube1Display_1.Assembly = ''
tube1Display_1.SelectOrientationVectors = 'Coordinates'
tube1Display_1.ScaleFactor = 0.1371263176202774
tube1Display_1.SelectScaleArray = 'Coordinates'
tube1Display_1.GlyphType = 'Arrow'
tube1Display_1.GlyphTableIndexArray = 'Coordinates'
tube1Display_1.GaussianRadius = 0.00685631588101387
tube1Display_1.SetScaleArray = ['POINTS', 'Coordinates']
tube1Display_1.ScaleTransferFunction = 'Piecewise Function'
tube1Display_1.OpacityArray = ['POINTS', 'Coordinates']
tube1Display_1.OpacityTransferFunction = 'Piecewise Function'
tube1Display_1.DataAxesGrid = 'Grid Axes Representation'
tube1Display_1.PolarAxes = 'Polar Axes Representation'
tube1Display_1.SelectInputVectors = ['POINTS', 'Coordinates']
tube1Display_1.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
tube1Display_1.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
tube1Display_1.ScaleTransferFunction.Points = [66.0, 0.0, 0.5, 0.0, 1981.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
tube1Display_1.OpacityTransferFunction.Points = [66.0, 0.0, 0.5, 0.0, 1981.0, 1.0, 0.5, 0.0]

# show data from threshold4
threshold4Display = Show(threshold4, renderView4, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold4Display.Representation = 'Surface'
threshold4Display.ColorArrayName = ['CELLS', 'PairType']
threshold4Display.LookupTable = pairTypeLUT
threshold4Display.PointSize = 50.0
threshold4Display.SelectTCoordArray = 'None'
threshold4Display.SelectNormalArray = 'None'
threshold4Display.SelectTangentArray = 'None'
threshold4Display.OSPRayScaleArray = 'Coordinates'
threshold4Display.OSPRayScaleFunction = 'Piecewise Function'
threshold4Display.Assembly = ''
threshold4Display.SelectOrientationVectors = 'Coordinates'
threshold4Display.ScaleFactor = 0.13713713139295577
threshold4Display.SelectScaleArray = 'Coordinates'
threshold4Display.GlyphType = 'Arrow'
threshold4Display.GlyphTableIndexArray = 'Coordinates'
threshold4Display.GaussianRadius = 0.006856856569647789
threshold4Display.SetScaleArray = ['POINTS', 'Coordinates']
threshold4Display.ScaleTransferFunction = 'Piecewise Function'
threshold4Display.OpacityArray = ['POINTS', 'Coordinates']
threshold4Display.OpacityTransferFunction = 'Piecewise Function'
threshold4Display.DataAxesGrid = 'Grid Axes Representation'
threshold4Display.PolarAxes = 'Polar Axes Representation'
threshold4Display.ScalarOpacityFunction = pairTypePWF
threshold4Display.ScalarOpacityUnitDistance = 0.013517097359597429
threshold4Display.OpacityArrayName = ['POINTS', 'Coordinates']
threshold4Display.SelectInputVectors = ['POINTS', 'Coordinates']
threshold4Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold4Display.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold4Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2047.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold4Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2047.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView5'
# ----------------------------------------------------------------

# show data from tube1
tube1Display_2 = Show(tube1, renderView5, 'GeometryRepresentation')

# trace defaults for the display properties.
tube1Display_2.Representation = 'Surface'
tube1Display_2.ColorArrayName = [None, '']
tube1Display_2.DiffuseColor = [0.0, 0.0, 0.0]
tube1Display_2.SelectTCoordArray = 'None'
tube1Display_2.SelectNormalArray = 'TubeNormals'
tube1Display_2.SelectTangentArray = 'None'
tube1Display_2.OSPRayScaleArray = 'Coordinates'
tube1Display_2.OSPRayScaleFunction = 'Piecewise Function'
tube1Display_2.Assembly = ''
tube1Display_2.SelectOrientationVectors = 'Coordinates'
tube1Display_2.ScaleFactor = 408.3247503757477
tube1Display_2.SelectScaleArray = 'Coordinates'
tube1Display_2.GlyphType = 'Arrow'
tube1Display_2.GlyphTableIndexArray = 'Coordinates'
tube1Display_2.GaussianRadius = 20.416237518787383
tube1Display_2.SetScaleArray = ['POINTS', 'Coordinates']
tube1Display_2.ScaleTransferFunction = 'Piecewise Function'
tube1Display_2.OpacityArray = ['POINTS', 'Coordinates']
tube1Display_2.OpacityTransferFunction = 'Piecewise Function'
tube1Display_2.DataAxesGrid = 'Grid Axes Representation'
tube1Display_2.PolarAxes = 'Polar Axes Representation'
tube1Display_2.SelectInputVectors = ['POINTS', 'Coordinates']
tube1Display_2.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
tube1Display_2.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
tube1Display_2.ScaleTransferFunction.Points = [0.00025549999554641545, 0.0, 0.5, 0.0, 301.7716064453125, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
tube1Display_2.OpacityTransferFunction.Points = [0.00025549999554641545, 0.0, 0.5, 0.0, 301.7716064453125, 1.0, 0.5, 0.0]

# show data from threshold5
threshold5Display = Show(threshold5, renderView5, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold5Display.Representation = 'Surface'
threshold5Display.ColorArrayName = ['CELLS', 'PairType']
threshold5Display.LookupTable = pairTypeLUT
threshold5Display.SelectTCoordArray = 'None'
threshold5Display.SelectNormalArray = 'None'
threshold5Display.SelectTangentArray = 'None'
threshold5Display.OSPRayScaleArray = 'Coordinates'
threshold5Display.OSPRayScaleFunction = 'Piecewise Function'
threshold5Display.Assembly = ''
threshold5Display.SelectOrientationVectors = 'Coordinates'
threshold5Display.ScaleFactor = 407.1
threshold5Display.SelectScaleArray = 'Coordinates'
threshold5Display.GlyphType = 'Arrow'
threshold5Display.GlyphTableIndexArray = 'Coordinates'
threshold5Display.GaussianRadius = 20.355
threshold5Display.SetScaleArray = ['POINTS', 'Coordinates']
threshold5Display.ScaleTransferFunction = 'Piecewise Function'
threshold5Display.OpacityArray = ['POINTS', 'Coordinates']
threshold5Display.OpacityTransferFunction = 'Piecewise Function'
threshold5Display.DataAxesGrid = 'Grid Axes Representation'
threshold5Display.PolarAxes = 'Polar Axes Representation'
threshold5Display.ScalarOpacityFunction = pairTypePWF
threshold5Display.ScalarOpacityUnitDistance = 114.96009771203582
threshold5Display.OpacityArrayName = ['POINTS', 'Coordinates']
threshold5Display.SelectInputVectors = ['POINTS', 'Coordinates']
threshold5Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
threshold5Display.OSPRayScaleFunction.Points = [37.35310363769531, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
threshold5Display.ScaleTransferFunction.Points = [4.023873805999756, 0.0, 0.5, 0.0, 506.97613525390625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
threshold5Display.OpacityTransferFunction.Points = [4.023873805999756, 0.0, 0.5, 0.0, 506.97613525390625, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = [renderView2, renderView4, renderView5]
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0

# ----------------------------------------------------------------
# restore active source
SetActiveSource(backpackpvd)

SaveScreenshot("ddmsExample.png", layout1)
