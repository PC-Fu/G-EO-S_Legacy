
from pylab import *
import glob
import os
from . import query, boundplot, colorplot, cxplot, clineplot, crplot, maxplot, sort_xdmv, decimate_xdmv


def process_all(path='*/', opt='2D'):
	pwd = os.getcwd()
	dirnames = sorted(glob.glob('%s/%s' % (pwd, path)))
	for d in dirnames:
		print('\n\n\nProcessing data in \'%s\' ...' % d)
		os.chdir(d)
		if (opt == '2D'):
			process_dir_2D()
		else:
			process_dir_3D()
	os.chdir(pwd)
	print '\n\nComplete!\n'


def process_dir_2D(export=True, res=True):
	if export:
		query.launch()
		query.zonal = ('FaceFields/Aperture')
		query.nodal = ('FaceNodeFields/nodalPressure')
		res = query.export_xdmv()
	if res:
		data = sort_xdmv()
		data = decimate_xdmv(data, dec=2)
		cxplot(data, 'FaceNodeFields/nodalPressure', xl=r'X ($m$)', yl=r'Pressure ($MPa$)', pt='', fname='pressure.pdf', scale=1.0e-6, showplot=False)
		maxplot(data, 'FaceNodeFields/nodalPressure', xl=r'T ($s$)', yl=r'Pressure ($MPa$)', pt='', fname='pmax.pdf', scale=1.0e-6, showplot=False)
		cxplot(data, 'FaceFields/Aperture', xl=r'X ($m$)', yl=r'Aperture ($\mu m$)', pt='', fname='aperture.pdf', scale=1.0e6, showplot=False)


def process_dir_3D(export=True, res=True):
	if export:
		query.launch()
		query.zonal = ('FaceFields/Aperture')
		query.nodal = ('FaceNodeFields/nodalPressure')
		res = query.export_xdmv()
	if res:
		data = sort_xdmv()
		crplot(data, 'FaceNodeFields/nodalPressure', xl=r'R ($m$)', yl=r'Pressure ($MPa$)', pt='', fname='pressure.pdf', scale=1.0e-6, showplot=False)
		clineplot(data, 'FaceNodeFields/nodalPressure', orient='y', xl=r'Y ($m$)', yl=r'Pressure ($MPa$)', pt='', fname='pressure2.pdf', scale=1.0e-6, showplot=False)
		maxplot(data, 'FaceNodeFields/nodalPressure', xl=r'T ($s$)', yl=r'Pressure ($MPa$)', pt='', fname='pmax.pdf', scale=1.0e-6, showplot=False)
		crplot(data, 'FaceFields/Aperture', xl=r'R ($m$)', yl=r'Aperture ($\mu m$)', pt='', fname='aperture.pdf', scale=1.0e6, showplot=False)
		clineplot(data, 'FaceFields/Aperture', orient='y', xl=r'Y ($m$)', yl=r'Aperture ($\mu m$)', pt='', fname='aperture2.pdf', scale=1.0e6, showplot=False)
		boundplot(data, orient=('z', 'y'), showplot=False)
		colorplot(data, 'FaceFields/Aperture', cycle=3221310, orient=('z', 'y'), cl=r'Aperture ($\mu m$)', pt='', fname='aperture_color.pdf', scale=1.0e6, showplot=False)
