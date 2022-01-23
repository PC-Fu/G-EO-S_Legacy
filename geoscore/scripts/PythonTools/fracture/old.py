
import glob
import os
from pylab import *
import sys
sys.path.append('/usr/gapps/visit/2.7.3/linux-x86_64/lib/site-packages/')

print 'Loading VisIt Python Interface...'
import visit as vs
vs.LaunchNowin()
vs.SuppressMessages(2)
vs.SuppressQueryOutputOn()
vs.OpenComputeEngine("localhost", ("-nn", "1", "-np", "1"))
print 'Ready!'


def sort_xdmv():
	# Read xdmv files into a list of dictionaries
	os.chdir('./xdmv')
	fnames = glob.glob('*.okc')
	fnames.sort()
	data = {}
	t = loadtxt('./time.csv', delimiter=',', unpack=True, skiprows=1)
	
	# Loop through files
	for f in fnames:
		name = f[:f.find('.')]

		# Build dictionary
		if name not in data.keys():
			xdmv = open(f, 'r')
			fsize = xdmv.readline()
			fsize.split()
			cycle = int(name[2:])
			index = where(t[0]==cycle)
			data[name] = {'keys': {}, 'cycle': cycle, 'time': t[1][index][0]}
			for ii in range(0, int(fsize[0])):
				key = xdmv.readline()[:-1]
				data[name]['keys'][ii] = key
				data[name][key] = array([])
			xdmv.close()

		# Load data
		tmp = loadtxt(f, delimiter='\t', unpack=True, skiprows=1+2*len(data[name]['keys']))
		for ii in range(0, len(tmp)):
			data[name][data[name]['keys'][ii]] = append(data[name][data[name]['keys'][ii]], array(tmp[ii]))

	# Sort results by first key (typically x)
	fnames = data.keys()
	for f in fnames:
		order = argsort(data[f][data[f]['keys'][0]])
		for k in data[f]['keys'].values():
			data[f][k] = data[f][k][order]

	# Separate nodal/zonal entries
	ndata = {}
	zdata = {}
	keys = sorted(data.keys())
	for k in keys:
		if (k[0] == 'n'):
			ndata[k[2:]] = data[k]
		else:
			zdata[k[2:]] = data[k]

	os.chdir('..')
	return ndata, zdata


def plot_all(ndata, zdata, showplot=True, dec=1):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	nkeys = sorted(ndata.keys())
	zkeys = sorted(zdata.keys())
	ignorekeys = ['cycle', 'time', 'keys', 'x', 'y', 'z']
	
	# Create a colormap
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=ndata[nkeys[0]]['time'], vmax=ndata[nkeys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Setup a plot for each appropriate key entry
	npkeys = sorted(ndata[nkeys[-1]].keys())
	zpkeys = sorted(ndata[zkeys[-1]].keys())
	pdict = {}
	ii = 0
	for p in (npkeys + zpkeys):
		if p not in ignorekeys:
			pdict[p] = {'handle': figure(ii), 'num': ii}
			ii += 1
			
	# Plot each entry as a function of distance
	for data in [ndata, zdata]:
		for k in data.keys():
			for p in data[k].keys():
				if p not in ignorekeys:
					figure(pdict[p].num)
					c = scalarMap.to_rgba(data[k]['time'])
					plot(data[k]['x'], data[k][p], color=c)
					hold(True)

	# Add annotations
	for p in pdict.keys():
		figure(pdict[p].num)
		xlabel(r'X (m)')
		ylabel(p)
		pdict[p].handle.add_axes([0.775, 0.45, 0.025, 0.4])
		cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
		cb.set_label('Time (s)')

	if showplot:
		show()
	else:
		for p in pdict.keys():
			pdict[p].handle.clf()


def plot_ap(ndata, zdata, showplot=True, dec=1):	
	from mpl_toolkits.axes_grid1 import host_subplot
	import mpl_toolkits.axisartist as AA
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	nkeys = sorted(ndata.keys())
	zkeys = sorted(zdata.keys())
	time = []
	aperture = []
	pressure = []

	# Create a colormap
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=ndata[nkeys[0]]['time'], vmax=ndata[nkeys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Aperture plot
	fig1 = figure(1)
	for ii in range(0, len(zkeys)):
		if (ii%dec == 0):
			k = zkeys[ii]
			time.append(zdata[k]['time'])
			aperture.append(max(zdata[k]['FaceFields/Aperture']))
			c = scalarMap.to_rgba(zdata[k]['time'])
			plot(zdata[k]['x'], 1e6*zdata[k]['FaceFields/Aperture'], color=c)
			hold(True)
	xlabel('X (m)')
	ylabel(r'Aperture ($\mu m$)')
	title('Fracture Propagation')
	axc = fig1.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/fracture_aperture.pdf', format='pdf')

	# Pressure plot
	fig2 = figure(2)
	for ii in range(0, len(nkeys)):
		if (ii%dec == 0):
			k = nkeys[ii]
			pressure.append(max(ndata[k]['FaceNodeFields/nodalPressure']))
			c = scalarMap.to_rgba(ndata[k]['time'])
			plot(ndata[k]['x'], 1e-6*ndata[k]['FaceNodeFields/nodalPressure'], color=c)
			hold(True)
	xlabel('X (m)')
	ylabel(r'Pressure ($MPa$)')
	title('Fracture Propagation')
	axc = fig2.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/fracture_pressure.pdf', format='pdf')

	# Max values plot
	fig3 = figure(3)
	host = host_subplot(111, axes_class=AA.Axes)
	subplots_adjust(left=0.175, right=0.9)
	par1 = host.twinx()
	p1, = host.plot(time, 1e6*array(aperture), 'b')
	p2, = par1.plot(time, 1e-6*array(pressure), 'r')
	host.set_xlabel(r'Time ($s$)')
	host.set_ylabel(r'Aperture ($\mu m$)')
	par1.set_ylabel(r'Pressure ($MPa$)')
	host.axis["left"].label.set_color(p1.get_color())
	par1.axis["right"].label.set_color(p2.get_color())
	title('Aperture vs. Pressure')
	legend([r'Max Aperture', r'Max Pressure'], loc=4)
	savefig('./fig/ap_time.pdf', format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()
		fig2.clf()
		fig3.clf()


def plot_max(ndata, zdata):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	nkeys = sorted(ndata.keys())
	zkeys = sorted(zdata.keys())
	time = []
	aperture = []
	pressure = []

	# Find max values
	for ii in range(0, len(zkeys)):
		k = zkeys[ii]
		time.append(zdata[k]['time'])
		aperture.append(max(zdata[k]['FaceFields/Aperture']))
		k = nkeys[ii]
		pressure.append(max(ndata[k]['FaceNodeFields/nodalPressure']))


	# Max values plot
	fig1 = figure(1)
	host = host_subplot(111, axes_class=AA.Axes)
	subplots_adjust(left=0.175, right=0.9)
	par1 = host.twinx()
	p1, = host.plot(time, 1e6*array(aperture), 'b')
	p2, = par1.plot(time, 1e-6*array(pressure), 'r')
	host.set_xlabel(r'Time ($s$)')
	host.set_ylabel(r'Aperture ($\mu m$)')
	par1.set_ylabel(r'Pressure ($MPa$)')
	host.axis["left"].label.set_color(p1.get_color())
	par1.axis["right"].label.set_color(p2.get_color())
	title('Aperture vs. Pressure')
	legend([r'Max Aperture', r'Max Pressure'], loc=4)
	savefig('./fig/ap_time.pdf', format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()


def export():
	if not os.path.isdir('./xdmv'):
		os.mkdir('./xdmv')

	# Set database
	pwd = os.getcwd()
	dbname = "localhost:%s/plot_* database" % (pwd)
	vs.OpenDatabase(dbname)

	# Setup fracture plot
	vs.AddPlot("Pseudocolor", "FaceFields/Aperture", 0, 0)
	vs.SetActivePlots(0)
	vs.AddOperator("Threshold", 0)
	ThresholdAtts = vs.ThresholdAttributes()
	ThresholdAtts.listedVarNames = ("FaceFields/Mass", "NodalFields/ghostRank")
	ThresholdAtts.lowerBounds = (1e-10, -1e+37)
	ThresholdAtts.upperBounds = (1e+37, -0.5)
	vs.SetOperatorOptions(ThresholdAtts, 0)
	vs.DrawPlots()
	vs.Query("SpatialExtents", use_actual_data=0)
	model_extents = vs.GetQueryOutputValue()
	flag_3D = ((model_extents[5]-model_extents[4]) > 0.0)

	ExportDBAtts = vs.ExportDBAttributes()
	ExportDBAtts.db_type = "Xmdv"
	ExportDBAtts.filename = "visit_ex_db"
	ExportDBAtts.dirname = "%s/" % pwd
	ExportDBAtts.variables = ("FaceFields/Aperture", "FaceFields/stressNOnFace")
	ExportDBAtts.opts.types = ()

	# Loop through database and export desired values
	n = vs.TimeSliderGetNStates()
	for ii in range(0, n):
		try:
			vs.SetTimeSliderState(ii)
			vs.Query("Cycle")
			cycle = vs.GetQueryOutputValue()
			vs.Query("Time")
			time = vs.GetQueryOutputValue()
			print 'VisIt: Cycle=%i, Time=%fs' % (cycle, time)
			ExportDBAtts.filename = "./xdmv/fracture_%i" % cycle
			vs.ExportDatabase(ExportDBAtts)
		except:
			pass

	vs.DeleteAllPlots()
	vs.CloseDatabase(dbname)


def wellbore():
	# Open database
	pwd = os.getcwd()
	dbname = "localhost:%s/plot_* database" % (pwd)
	vs.OpenDatabase(dbname)

	# Setup wellbore plot
	vs.AddPlot("Pseudocolor", 'FaceNodeFields/nsigma_x', 0, 0)
	vs.SetActivePlots(0)
	vs.DrawPlots()
	vs.Query("SpatialExtents", use_actual_data=0)
	model_extents = vs.GetQueryOutputValue()

	# Add threshold/box operator
	vs.AddOperator("Threshold", 0)
	ThresholdAtts = vs.ThresholdAttributes()
	ThresholdAtts.listedVarNames = ("FaceFields/isExternal", "NodalFields/ghostRank", "FaceFields/ruptureState")
	ThresholdAtts.lowerBounds = (0.5, -1e+37, -1e+37)
	ThresholdAtts.upperBounds = (1e+37, -0.5, 0.5)
	vs.AddOperator("Box", 0)
	BoxAtts = vs.BoxAttributes()
	BoxAtts.minx = 0.99*model_extents[0]
	BoxAtts.maxx = 0.99*model_extents[1]
	BoxAtts.miny = 0.99*model_extents[2]
	BoxAtts.maxy = 0.99*model_extents[3]
	BoxAtts.minz = 0.99*model_extents[4]
	BoxAtts.maxz = 0.99*model_extents[5]
	vs.SetOperatorOptions(BoxAtts, 0)
	vs.DrawPlots()

	# Find target zone for tracking wellbore pressure
	vs.Query("SpatialExtents", use_actual_data=1)
	extents = vs.GetQueryOutputValue()
	wellbore_target = (extents[0], 0.5*(extents[2]+extents[3]), extents[5])
	vs.ZonePick(coord=(extents[0], 0.5*(extents[2]+extents[3]), extents[5]))
	ptmp = vs.GetPickOutputObject()
	wellbore_zone = ptmp['zone_id']
	if ((model_extents[5]-model_extents[4]) > 0.0):
		wellbore_domain = ptmp['domain_id']
	else:
		wellbore_domain = 1

	# Setup output file
	f = open('wellbore_pressure.csv', 'w')
	f.write('cycle, time (s), pressure (Pa)\n')

	# Query over time
	n = vs.TimeSliderGetNStates()
	print '\n\nGathering wellbore pressure data from %i plots:' % n
	bar = progress.ProgressBar()
	for ii in range(0, n):
		bar.update(float(ii)/n)
		try:
			# Advance slider
			vs.SetActivePlots(0)
			vs.SetTimeSliderState(ii)
			vs.Query("Cycle")
			cycle = vs.GetQueryOutputValue()
			vs.Query("Time")
			time = vs.GetQueryOutputValue()
			
			# Find pressure
			vs.PickByZone(wellbore_zone, wellbore_domain)
			ptmp = vs.GetPickOutputObject()
			pressure = -1*min(ptmp['FaceNodeFields/nsigma_x'].values())
			f.write('%i, %f, %1.4e\n' %(cycle, time, pressure))
		except:
			pass

	bar.update(1.0)
	f.close()
	vs.DeleteAllPlots()
	vs.CloseDatabase(dbname)


def query():
	# Set database
	pwd = os.getcwd()
	dbname = "localhost:%s/plot_* database" % (pwd)
	vs.OpenDatabase(dbname)

	# Setup fracture plot
	vs.AddPlot("Pseudocolor", "FaceFields/Aperture", 0, 0)
	vs.SetActivePlots(0)
	vs.AddOperator("Threshold", 0)
	ThresholdAtts = vs.ThresholdAttributes()
	ThresholdAtts.outputMeshType = 0
	ThresholdAtts.listedVarNames = ("FaceFields/ruptureState", "NodalFields/ghostRank")
	ThresholdAtts.zonePortions = (1, 1)
	ThresholdAtts.lowerBounds = (1.5, -1e+37)
	ThresholdAtts.upperBounds = (1e+37, -0.5)
	ThresholdAtts.defaultVarName = "FaceFields/Aperture"
	ThresholdAtts.defaultVarIsScalar = 0
	vs.SetOperatorOptions(ThresholdAtts, 0)
	vs.DrawPlots()
	vs.Query("SpatialExtents", use_actual_data=0)
	model_extents = vs.GetQueryOutputValue()
	flag_3D = ((model_extents[5]-model_extents[4]) > 0.0)


	# Setup output csv file
	f = open('spatial_query.csv', 'w')
	f.write('cycle, time, fx, fy, fz, xmin, xmax, ymin, ymax, zmin, zmax, area, volume, mass, amax, aavg, p_wb, p_fr, qmax\n')

	# Go to starting time
	n = vs.TimeSliderGetNStates()
	for ii in range(0, n):
		try:
			# Geometric quantities
			vs.SetActivePlots(0)
			vs.SetTimeSliderState(ii)
			vs.Query("Cycle")
			cycle = vs.GetQueryOutputValue()
			vs.Query("Time")
			time = vs.GetQueryOutputValue()
			print 'VisIt: Cycle=%i, Time=%fs' % (cycle, time)
			vs.Query("SpatialExtents", use_actual_data=1)
			extents = vs.GetQueryOutputValue()
			vs.Query("3D surface area")
			area = vs.GetQueryOutputValue()
			vs.ChangeActivePlotsVar("FaceFields/Volume")
			vs.Query("Variable Sum", use_actual_data=1)
			volume = vs.GetQueryOutputValue()

			# Fracture tip
			if flag_3D:
				vs.ChangeActivePlotsVar("FaceFields/birthTime")
			else:
				vs.ChangeActivePlotsVar("NodalFields/SIF_II")
			vs.Query("Max", use_actual_data=0)
			fractureTip = vs.GetQueryOutputObject()

			# Mass, Aperature, pressure
			vs.Query("NumNodes", use_actual_data=1)
			N = vs.GetQueryOutputValue()
			vs.Query("NumZones", use_actual_data=1)
			M = vs.GetQueryOutputValue()
			vs.ChangeActivePlotsVar("FaceFields/Mass")
			vs.Query("Variable Sum", use_actual_data=1)
			mass = vs.GetQueryOutputValue()
			
			#Alternate pressure estimate
			vs.ChangeActivePlotsVar("FaceFields/stressNOnFace")
			vs.Query("Min", use_actual_data=1)
			pressure = -vs.GetQueryOutputValue()
			vs.Query("Variable Sum", use_actual_data=1)
			pressure_avg = -vs.GetQueryOutputValue() / N
			# vs.ChangeActivePlotsVar("FaceNodeFields/nodalPressure")
			# vs.Query("Max", use_actual_data=1)
			# pressure = vs.GetQueryOutputValue()
			# vs.Query("Variable Sum", use_actual_data=1)
			# pressure_avg = vs.GetQueryOutputValue() / N
			vs.ChangeActivePlotsVar("FaceFields/flowRate")
			vs.Query("Max", use_actual_data=1)
			qmax = vs.GetQueryOutputValue()
			vs.ChangeActivePlotsVar("FaceFields/Aperture")
			vs.Query("Max", use_actual_data=1)
			aperture = vs.GetQueryOutputValue()
			vs.ChangeActivePlotsVar("FaceFields/Aperture")
			vs.Query("Variable Sum", use_actual_data=1)
			aperture_avg = vs.GetQueryOutputValue() / M

			# Write data to csv file
			f.write('%i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %1.4e, %1.4e, %1.4e, %1.4e, %1.4e, %1.4e, %1.4e\n' 
					%(cycle, time, fractureTip['max_coord'][0], fractureTip['max_coord'][1], fractureTip['max_coord'][2],
					  extents[0], extents[1], extents[2], extents[3], extents[4], extents[5], 
				      area, volume, mass, aperture, aperture_avg, pressure, pressure_avg, qmax))
		except:
			print 'Query Error: Cycle=%i, Time=%fs' % (cycle, time)

	f.close()
	vs.DeleteAllPlots()
	vs.CloseDatabase(dbname)
	



def plot_query(fname='spatial_query.csv', ptitle='Hydrofracture Simulation', showplot=True):
	from mpl_toolkits.axes_grid1 import host_subplot
	import mpl_toolkits.axisartist as AA

	# Read data
	cycle, t, fx, fy, fz, xmin, xmax, ymin, ymax, zmin, zmax, area, volume, mass, amax, aavg, pmax, pavg, qmax = loadtxt(fname, delimiter=',', unpack=True, skiprows=1)
	r_o = 1
	L = sqrt(area/pi + r_o**2) - r_o
	pmax *= 1e-6
	pavg *= 1e-6
	amax *= 1e6
	aavg *= 1e6

	# Plot extents
	figure()
	plot(t, xmin, 'b')
	hold(True)
	plot(t, xmax, 'b--')
	plot(t, ymin, 'r')
	plot(t, ymax, 'r--')
	plot(t, zmin, 'g')
	plot(t, zmax, 'g--')
	xlabel(r'Time ($s$)')
	ylabel(r'Distance ($m$)')
	title(ptitle + '\nExtents')
	legend([r'$x_{min}$', r'$x_{max}$', r'$y_{min}$', r'$y_{max}$', r'$z_{min}$', r'$z_{max}$', ], loc=1)
	savefig('fracture_extents.pdf', format='pdf')

	# Plot aperture, pressure vs. time
	figure()
	host = host_subplot(111, axes_class=AA.Axes)
	subplots_adjust(left=0.175, right=0.9)
	par1 = host.twinx()
	p1, = host.plot(t, amax, 'b')
	p2, = par1.plot(t, pmax, 'r')
	par1.hold(True)
	p3, = par1.plot(t, pavg, 'r--')
	host.set_xlabel(r'Time ($s$)')
	host.set_ylabel(r'Aperture ($\mu m$)')
	par1.set_ylabel(r'Pressure ($MPa$)')
	host.axis["left"].label.set_color(p1.get_color())
	par1.axis["right"].label.set_color(p2.get_color())
	title(ptitle + '\nAperture vs. Pressure')
	legend([r'Aperture', r'Wellbore Pressure', r'Fracture Pressure'], loc=4)
	savefig('ap_time.pdf', format='pdf')


	# Fracture Trace
	# Filter out bad points
	fx = insert(fx, 0, 0.02)
	fy = insert(fy, 0, 0.1)
	filt = []
	for ii in range(0, len(fx)):
		if ((fx[ii]>=fx[0]) and (fy[ii]>=fy[0])):
			filt.append(ii)

	figure()
	tmp = linspace(0, 2*pi, 1000)
	plot(0.1*sin(tmp), 0.1*cos(tmp), 'k')
	hold(True)
	plot(fx[filt], fy[filt], 'b')
	xlabel(r'X ($m$)')
	ylabel(r'Y ($m$)')
	xlim([-1.0, 1.0])
	ylim([-1.0, 1.0])
	title(ptitle + '\nFracture Trace')
	axis('equal')
	legend([r'Borehole', r'Fracture Path'], loc=4)
	savefig('fracture_trace.pdf', format='pdf')

	if showplot:
		show()



def plot_multiple(showplot=True):
	# Read data from *csv files into a list of dicts.
	fnames = glob.glob('*.csv')
	fnames.sort()
	data = []
	for f in fnames:
		cycle, t, fx, fy, fz, xmin, xmax, ymin, ymax, zmin, zmax, area, volume, mass, amax, aavg, pmax, pavg, qmax = loadtxt(f, delimiter=',', unpack=True, skiprows=1)
		data.append({'name': f[:-4],
					 'cycle': cycle,
					 'time': t,
					 'path': [fx, fy, fz],
					 'limits': [xmin, xmax, ymin, ymax, zmin, zmax],
					 'geometric': [area, volume, mass],
					 'aperture': amax*1e6,
					 'pressure': pmax*1e-6,
					 'flow': qmax,
					 'filter': [],
					 'offset': 0.0})

	style = ['b', 'c', 'g', 'y', 'r', 'm', 
			 'b--', 'c--', 'g--', 'y--', 'r--', 'm--',
			 'b:', 'c:', 'g:', 'y:', 'r:', 'm:']

	#style = ['c', 'c--', 'b', 'r']

	# # Estimated values for chi=0.35 (Zhang et al. 2011)
	f2t = array([0.0, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3])/10
	f2p = array([0.0, 38.0, 38.0, 34.5, 32.0, 31.0, 30.4, 29.2, 28.7, 28.0])
	f2x = array([0.02, 0.035, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30])
	f2y = array([0.1, 0.14, 0.19, 0.218, 0.235, 0.245, 0.255, 0.261, 0.265])


	# Filter out incorrect fracture path points, find offset value
	pref = 10.0
	px_range = [0.01, 10]
	py_range = [0.0, 10]
	for ii in range(0, len(data)):
		pind = find(data[ii]['pressure'] < pref)
		data[ii]['offset'] = data[ii]['time'][pind[-1]]

		for jj in range(0, len(data[ii]['time'])):
			if ((px_range[0] <= data[ii]['path'][0][jj] <= px_range[1]) &  (py_range[0] <= data[ii]['path'][1][jj] <= py_range[1])):
				data[ii]['filter'].append(jj)

	# Fracture path
	figure()
	tmp = linspace(0, 2*pi, 1000)
	plot(0.1*sin(tmp), 0.1*cos(tmp), 'k', label='wellbore')
	hold(True)
	for ii in range(0, len(data)):
		f = data[ii]['filter']
		plot(data[ii]['path'][0][f], data[ii]['path'][1][f], style[ii], label=data[ii]['name'])
	plot(f2x, f2y, style[ii+1], label='Zhang et al. (2011)')
	xlim([-.1, 0.4])
	ylim([-.1, 0.4])
	xlabel(r'X ($m$)')
	ylabel(r'Y ($m$)')
	legend(loc=4)
	title('Fracture Path')
	savefig('fracture_path.pdf', format='pdf')

	# Wellbore pressure
	figure()
	hold(True)
	for ii in range(0, len(data)):
		plot(data[ii]['time'] - data[ii]['offset'], data[ii]['pressure'], style[ii], label=data[ii]['name'])
	plot(f2t, f2p, style[ii+1], label='Zhang et al. (2011)')
	xlabel(r'Time ($s$)')
	ylabel(r'Pressure ($MPa$)')
	legend(loc=4)
	title('Wellbore Pressure vs. Time')
	savefig('wellbore_pressure.pdf', format='pdf')

	# Fracture aperture
	figure()
	hold(True)
	for ii in range(0, len(data)):
		plot(data[ii]['time'] - data[ii]['offset'], data[ii]['aperture'], style[ii], label=data[ii]['name'])
	xlabel(r'Time ($s$)')
	ylabel(r'Aperture ($\mu m$)')
	legend(loc=4)
	title('Wellbore Aperture vs. Time')
	savefig('wellbore_aperture.pdf', format='pdf')

	if showplot:
		show()


