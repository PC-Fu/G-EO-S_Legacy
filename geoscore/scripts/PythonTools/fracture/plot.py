
import os
from pylab import *
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar


def colorplot(data, pkey, cycle=0, orient=('x', 'y'), cl=r'C', pt='', fname='plot.pdf', scale=1, showplot=True):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	keys = sorted(data.keys())

	

	# Plot values vs. X
	fig1 = figure(1)
	for k in keys:
		if pkey in data[k].keys():
			if (data[k]['cycle'] == cycle):
				xi = linspace(min(data[k][orient[0]]), max(data[k][orient[0]]), 200)
				yi = linspace(min(data[k][orient[1]]), max(data[k][orient[1]]), 200)
				zi = griddata(data[k][orient[0]], data[k][orient[1]], data[k][pkey], xi, yi, interp='linear')
				pcolormesh(xi, yi, zi*scale)
				cb = colorbar()
				cb.set_label(cl)
				# cs = contour(xi, yi, zi*scale, 10)
				# clabel(cs, inline=1, fontsize=10, fmt=r'%i $\mu m$')

	# Add annotations and save
	xlabel(r'%s ($m$)' % orient[0])
	ylabel(r'%s ($m$)' % orient[1])
	title(pt)
	axis('equal')
	savefig('./fig/%s' % fname, format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()


def cxplot(data, pkey, xl=r'X ($m$)', yl=r'Y', pt='', fname='plot.pdf', scale=1, showplot=True):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	keys = sorted(data.keys())

	# Build colormap
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=data[keys[0]]['time'], vmax=data[keys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Plot values vs. X
	fig1 = figure(1)
	for k in keys:
		for p in data[k].keys():
			if (p == pkey):
				c = scalarMap.to_rgba(data[k]['time'])
				plot(data[k]['x'], data[k][p]*scale, color=c)
				hold(True)

	# Add annotations and save
	xlabel(xl)
	ylabel(yl)
	title(pt)
	axc = fig1.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/%s' % fname, format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()


def crplot(data, pkey, xl=r'R ($m$)', yl=r'Y', pt='', fname='plot.pdf', scale=1, showplot=True):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	keys = sorted(data.keys())

	# Build colormap
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=data[keys[0]]['time'], vmax=data[keys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Plot values vs. radius
	fig1 = figure(1)
	for k in keys:
		for p in data[k].keys():
			if (p == pkey):
				c = scalarMap.to_rgba(data[k]['time'])
				r = array(sqrt(data[k]['x']**2 + data[k]['y']**2 + data[k]['z']**2))
				d = array(data[k][p]*scale)
				order = argsort(r)
				r = r[order]
				d = d[order]
				
				# Apply moving max filter
				dr = max([1.1*max(diff(r)), 0.001])
				r2, d2 = mov_max(r, d, dr)

				plot(r2, d2, color=c)
				hold(True)

	# Add annotations and save
	xlabel(xl)
	ylabel(yl)
	title(pt)
	axc = fig1.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/%s' % fname, format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()



def clineplot(data, pkey, orient='x', xl=r'L ($m$)', yl=r'Y', pt='', fname='plot.pdf', scale=1, showplot=True):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	keys = sorted(data.keys())

	# Build colormap
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=data[keys[0]]['time'], vmax=data[keys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Plot values vs. radius
	fig1 = figure(1)
	for k in keys:
		for p in data[k].keys():
			if (p == pkey):
				c = scalarMap.to_rgba(data[k]['time'])
				r = array(data[k][orient])
				d = array(data[k][p]*scale)
				order = argsort(r)
				r = r[order]
				d = d[order]
				
				# Apply moving max filter
				dr = max([1.1*max(diff(r)), 0.001])
				r2, d2 = mov_max(r, d, dr)

				plot(r2, d2, color=c)
				hold(True)

	# Add annotations and save
	xlabel(xl)
	ylabel(yl)
	title(pt)
	axc = fig1.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/%s' % fname, format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()


def boundplot(data, orient=('x', 'y'), showplot=True):
	# Check version of scipy
	import scipy
	if (scipy.version.version < '0.12.0'):
		warning.warn('Scipy v0.12.0 required (on LC use python 2.7.7 bin)')

	from scipy.spatial import ConvexHull
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')

	# Create a colormap
	keys = sorted(data.keys())
	jet = cm.get_cmap('jet')
	cNorm = mcolors.Normalize(vmin=data[keys[0]]['time'], vmax=data[keys[-1]]['time'])
	scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

	# Boundary plot
	fig1 = figure(1)
	for k in keys:
		if (k[0] == 'n'):
			c = scalarMap.to_rgba(data[k]['time'])
			frac = swapaxes(vstack((data[k][orient[0]], data[k][orient[1]])), 0, 1)
			hull = ConvexHull(frac)
			fill(frac[hull.vertices, 0], frac[hull.vertices, 1], edgecolor=c, fill=False)
			hold(True)
	xlabel(r'%s ($m$)' % orient[0])
	ylabel(r'%s ($m$)' % orient[1])
	title('Fracture Extents')
	axis('equal')
	# ylim([0, 1])
	axc = fig1.add_axes([0.775, 0.45, 0.025, 0.4])
	cb = mcolorbar.ColorbarBase(axc, norm=cNorm, cmap=jet)
	cb.set_label('Time (s)')
	savefig('./fig/extents.pdf', format='pdf')
	
	if showplot:
		show()
	else:
		fig1.clf()


def maxplot(data, pkey, xl=r'Time ($s$)', yl=r'Y', pt='', fname='plot.pdf', scale=1, showplot=True):
	# Setup
	if not os.path.isdir('./fig'):
		os.mkdir('./fig')
	keys = sorted(data.keys())

	# Collect t, pkey values
	t = []
	g = []
	for k in keys:
		for p in data[k].keys():
			if (p == pkey):
				t.append(data[k]['time'])
				g.append(max(data[k][p])*scale)

	fig1 = figure(1)
	plot(t, g, 'b')
	xlabel(xl)
	ylabel(yl)
	title(pt)
	savefig('./fig/%s' % fname, format='pdf')

	if showplot:
		show()
	else:
		fig1.clf()


def mov_max(x, y, dx):
	# Setup
	n = int((max(x) - min(x))/dx)
	x2 = linspace(min(x), max(x), n)
	y2 = array(zeros((n)))

	for ii in range(0, n):
		# Find max value in window
		ind = (x2[ii]-0.5*dx <= x) & (x <= x2[ii]+0.5*dx)
		y2[ii] = max(y[ind])
		
		# Adjust x position
		tmp = where((y==y2[ii]) & (x2[ii]-0.5*dx <= x) & (x <= x2[ii]+0.5*dx))
		x2[ii] = x[tmp[0][0]] 

	return x2, y2