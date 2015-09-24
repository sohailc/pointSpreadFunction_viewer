"""
This code calculates the point spread function (PSF) of an wide field optical
microscope and renders the resulting image in 3D via the volumetric 
rendering function in VTK. 

The point spread function of an optical system is the intensity field 
which is generated when an ideal point source of light is imaged. 

The function which calculates the PSF is implemented in C and is accessed
in Python via ctypes. 

"""


import numpy
import ctypes
import cPickle
import os
from vtkVolumeRender import main as volumeRender

class microscopeParams(ctypes.Structure):
    _fields_ = [("NA", ctypes.c_double), \
                ("NoilReal", ctypes.c_double), \
                ("NoilDesign", ctypes.c_double),\
                ("NcovReal", ctypes.c_double),\
                ("NcovDesign", ctypes.c_double),\
                ("NspecReal", ctypes.c_double),\
                ("NspecDesign", ctypes.c_double),\
                ("ToilReal", ctypes.c_double),\
                ("ToilDesign", ctypes.c_double),\
                ("TcovReal", ctypes.c_double),\
                ("TcovDesign", ctypes.c_double), \
                ("TspecReal", ctypes.c_double),\
                ("TspecDesign", ctypes.c_double),\
                ("zdReal", ctypes.c_double), \
                ("zdDesign", ctypes.c_double), \
                ("M", ctypes.c_double)]

def convertToUInt8(data):
	
	mx = numpy.max(numpy.ravel(data))
	mn = numpy.min(numpy.ravel(data))
	
	# map data: min = 0, max = 255
	
	n = 255*(data - mn) / (mx - mn)
	
	return n.astype(numpy.uint8)

def generatePSF():

	# we setup the properties of the microscope. 
	params = {	"NA": 			1.4,		# numerical aparture of the lens
				"NoilReal": 	1.515,		# the actual refrective index of the oil
				"NoilDesign":	1.515,		# the nominal refractive index of the oil
				"NcovReal":		1.530,		# the actual refractive index if the cover slip
				"NcovDesign":	1.530,		# the nominal refractive index if the cover slip
				"NspecReal":	1.33,		# the actual refrective index of the specimin
				"NspecDesign":	1.33,		# the nominal refrective index of the specimin
				"ToilReal":		0.25e-3,	# the actual thickness of the oil
				"ToilDesign":	0.25e-3,	# the nominal thickness of the oil
				"TcovReal":		170e-6,		# the actual thickness of the cover slip
				"TcovDesign":	170e-6,		# the nominal thickness of the cover slip
				"TspecReal":	3.0E-6,		# the actual thickness of the specimin 
				"TspecDesign":	0.0E-6,		# the nominal thickness of the specimin 
				"zdReal":		16e-2,		# the actual distance between the projection plane and the apature 
				"zdDesign":		16e-2,		# the nominal distance between the projection plane and the apature
				"M":  			100}		# the magnification of the lens
	
	paramsC = microscopeParams(**params)
	
	sizeZ = 2*128
	sizeXY = 2*64
	wavelength = ctypes.c_double(500E-9)
	oversamplingFactor = ctypes.c_double(4.0)
	imageC = ctypes.POINTER(ctypes.c_double)()
	scales = ctypes.POINTER(ctypes.c_double)()
	
	psfLib = ctypes.CDLL("./confocalPSF.so")
	calculate3dPSF = psfLib.calculate3dPSF
	
	print "calculating point spread function"
	calculate3dPSF(	paramsC, 
					sizeZ, 
					sizeXY, 
					wavelength, 
					oversamplingFactor, 
					ctypes.byref(imageC), 
					ctypes.byref(scales))
	
	imSize = sizeZ*sizeXY**2
	
	image = numpy.array(imageC[:imSize])
	image = numpy.reshape(image, (sizeZ,sizeXY,sizeXY))
	
	image = convertToUInt8(image) # only 8 bit unsigned integer 
	# data can be shown in VTK via volume rendering. 
	
	return image, params, wavelength.value
	
def main():
	
	if not os.path.exists("PSF.pickle"):
		image, params, wavelength = generatePSF()
		
		with open("PSF.pickle","wb") as fh:
			cPickle.dump((image, params, wavelength), fh)
		
	else:
		with open("PSF.pickle","rb") as fh:
			image, params, wavelength = cPickle.load(fh)
	
	description = "Point Spread Function simulation with lambda = %.1f nm" % (wavelength*1E9)
	
	volumeRender(image, description)

if __name__ == "__main__":
	main()
	
	
    
	
	
