"""
This code has been adapted from:
http://www.vtk.org/Wiki/VTK/Examples/Python/vtkWithNumpy
"""

import vtk
import numpy as np

def alphaChannelFunction():
	f = vtk.vtkPiecewiseFunction()
	f.AddPoint(0, 0.0)
	f.AddPoint(50, 0.05)
	f.AddPoint(100, 0.2)
	f.AddPoint(150, 0.6)
	f.AddPoint(255, 0.7)
	return f

def colorChannelFunction():
	
	f = vtk.vtkColorTransferFunction()
	f.AddRGBPoint(0, 0.0, 0.0, 1.0)
	f.AddRGBPoint(255, 1.0, 1.0, 1.0)
	
	return f

def createVolumeProperty():
	
	alphaChannelFunc = alphaChannelFunction()
	colorChannelFunc = colorChannelFunction()
	
	volumeProperty = vtk.vtkVolumeProperty()
	volumeProperty.SetColor(colorChannelFunc)
	volumeProperty.SetScalarOpacity(alphaChannelFunc)
	
	return volumeProperty
	
def createVolumeMapper():
	
	compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
	volumeMapper = vtk.vtkVolumeRayCastMapper()
	volumeMapper.SetVolumeRayCastFunction(compositeFunction)
	
	return volumeMapper

def createRenderWindow(volume, backgroundColor, size):
	
	renderer = vtk.vtkRenderer()

	renderer.AddVolume(volume)
	renderer.SetBackground(*backgroundColor)
	
	renderWin = vtk.vtkRenderWindow()
	renderWin.AddRenderer(renderer)
	renderWin.SetSize(*size)
	renderInteractor = vtk.vtkRenderWindowInteractor()
	renderInteractor.SetRenderWindow(renderWin)
	
	return renderer, renderWin

def addText(renderer, text):
	textActor = vtk.vtkTextActor()
	textActor.GetTextProperty().SetFontSize (24)
	textActor.SetPosition2(10, 40)
	renderer.AddActor2D(textActor)
	textActor.SetInput (text)
	textActor.GetTextProperty().SetColor( 1.0,0.0,0.0 )

def addAxis(renderWindowInteractor):
	
	axis = vtk.vtkAxesActor()
	widget = vtk.vtkOrientationMarkerWidget()
	widget.SetOrientationMarker(axis)
	widget.SetInteractor(renderWindowInteractor)
	widget.EnabledOn()
	widget.InteractiveOn()
	

def exitCheck(obj, event):
    if obj.GetEventPending() != 0:
        obj.SetAbortRender(1)

def main(data, description=""): # data is a nD numpy array 
	
	m,l,k = data.shape
	
	dataImporter = vtk.vtkImageImport()
	dataString = data.tostring()
	dataImporter.CopyImportVoidPointer(dataString, len(dataString))
	dataImporter.SetDataScalarTypeToUnsignedChar()
	dataImporter.SetNumberOfScalarComponents(1)
	
	dataImporter.SetDataExtent(0, k-1, 0, l-1, 0, m-1)
	dataImporter.SetWholeExtent(0, k-1, 0, l-1, 0, m-1)
	
	volumeProperty = createVolumeProperty()
	volumeMapper = createVolumeMapper()
	volumeMapper.SetInputConnection(dataImporter.GetOutputPort())
	
	volume = vtk.vtkVolume()
	volume.SetMapper(volumeMapper)
	volume.SetProperty(volumeProperty)
	volume.SetOrientation(90,0,0)
	
	renderer, renderWin = createRenderWindow(volume, (0,0,0), (800, 800))
	addText(renderer, description)
	
	renderWin.AddObserver("AbortCheckEvent", exitCheck)
	
	renderInteractor = vtk.vtkRenderWindowInteractor()
	renderInteractor.SetRenderWindow(renderWin)
	#addAxis(renderInteractor)
	#axis = vtk.vtkAxesActor()
	#renderer.AddActor(axis)
	renderer.ResetCamera()
 
	renderInteractor.Initialize()

	renderWin.Render()
	renderInteractor.Start()
	
