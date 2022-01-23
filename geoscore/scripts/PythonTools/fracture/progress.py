#!/usr/local/bin/python
import sys

class ProgressBar:
	def __init__(self):
		sys.stdout.write('\n[%-20s] %d%%' % ('='*0, 5*0))
		sys.stdout.flush()

	def update(self, percent):
		ii = int(20*percent)
		sys.stdout.write('\r')
		sys.stdout.write('[%-20s] %d%%' % ('='*ii, 5*ii))
		sys.stdout.flush()
