
import os
from pylab import *
import sys
from . import ProgressBar

class Query():
	def __init__(self):
		self.engine = False
		self.dbname = 'plot_'
		self.zonal = ("FaceFields/Aperture")
		self.nodal = ("FaceNodeFields/nodalPressure")

	def launch(self):
		if not self.engine:
			sys.path.append('/usr/gapps/visit/2.8.1/linux-x86_64/lib/site-packages/')
			print '\n\nLoading VisIt Python Interface...\n'
			self.vs = __import__('visit')
			self.vs.LaunchNowin()
			self.vs.SuppressMessages(1)
			self.vs.SuppressQueryOutputOn()
			# self.vs.OpenComputeEngine("localhost", ("-nn", "1", "-np", "1"))
			print '\nReady!\n'

			self.engine = True

	def export_xdmv(self):
		# Set database
		self.launch()
		pwd = os.getcwd()
		dbname = "localhost:%s/%s* database" % (pwd, self.dbname)
		res = self.vs.OpenDatabase(dbname)

		if res:
			# Setup fracture plot
			self.vs.AddPlot("Pseudocolor", "FaceFields/Aperture", 0, 0)
			self.vs.SetActivePlots(0)
			self.vs.AddOperator("Threshold", 0)
			ThresholdAtts = self.vs.ThresholdAttributes()
			ThresholdAtts.listedVarNames = ("FaceFields/flowFaceType", "NodalFields/ghostRank")
			ThresholdAtts.lowerBounds = (-0.5, -1e+37)
			ThresholdAtts.upperBounds = (0.0, -0.5)
			self.vs.SetOperatorOptions(ThresholdAtts, 0)
			self.vs.DrawPlots()

			ExportDBAtts = self.vs.ExportDBAttributes()
			ExportDBAtts.db_type = "Xmdv"
			ExportDBAtts.filename = "visit_ex_db"
			ExportDBAtts.dirname = "%s/" % pwd
			ExportDBAtts.variables = ("FaceFields/Aperture")
			ExportDBAtts.opts.types = ()

			# Record cycle/time information
			if not os.path.isdir('./xdmv'):
				os.mkdir('./xdmv')
			if os.path.isfile('./xdmv/time.csv'):
				f = open('./xdmv/time.csv', 'a')
			else:
				f = open('./xdmv/time.csv', 'a')
				f.write('cycle, time\n')

			# Loop through database and export desired values
			n = self.vs.TimeSliderGetNStates()
			print 'Gathering fracture data from %i plots:' % n
			bar = ProgressBar()
			for ii in range(0, n):
				try:
					bar.update(float(ii)/n)
					self.vs.SetTimeSliderState(ii)
					self.vs.Query("Cycle")
					cycle = self.vs.GetQueryOutputValue()
					self.vs.Query("Time")
					time = self.vs.GetQueryOutputValue()
					f.write('%i, %f\n' % (cycle, time))

					if self.zonal:
						ExportDBAtts.variables = self.zonal
						ExportDBAtts.filename = "./xdmv/z_%08d" % cycle
						self.vs.ExportDatabase(ExportDBAtts)

					if self.nodal:
						ExportDBAtts.variables = self.nodal
						ExportDBAtts.filename = "./xdmv/n_%08d" % cycle
						self.vs.ExportDatabase(ExportDBAtts)
				except:
					pass
			
			bar.update(1.0)
			self.vs.DeleteAllPlots()
			self.vs.CloseDatabase(dbname)
			f.close()

		return res

query = Query()