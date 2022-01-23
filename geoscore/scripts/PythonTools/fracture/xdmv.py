
from pylab import *
import os
import glob


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

	os.chdir('..')
	return data


def decimate_xdmv(data, dec=2):
	data2 = data.copy()

	# Separate nodal and zonal keys
	keys = sorted(data2.keys())
	nkeys = []
	zkeys = []
	for k in keys:
		if (k[0] == 'n'):
			nkeys.append(k)
		else:
			zkeys.append(k)

	for ii in range(0, len(nkeys)):
		if (ii%dec != 0):
			del data2[nkeys[ii]]

	for ii in range(0, len(zkeys)):
		if (ii%dec != 0):
			del data2[zkeys[ii]]

	return data2