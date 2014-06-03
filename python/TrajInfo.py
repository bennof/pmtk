'''

    ./python/TrajInfo.py 

    Protein Motion TK (pmtk) - A library and several tools to estimate protein motions from varying input
    Copyright (C) 2009-2014  Benjamin Falkner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


'''



import Tkinter as tk
from TkTools import *


import Plot





class TrajInfo():
	def __init__(self,parent):
		self.view=TrajInfoView(self,parent)
		pass
	def onPlot(self):
		pass
	def onUpdate(self,data):
		self.view.update(data)
	def onPopUp(self,e=None):
		print " (%i,%i)"%(e.x_root, e.y_root)
		try:
			self.view.popupm.tk_popup(e.x_root, e.y_root)
		except Exception as err:
			exc_type, exc_value, exc_traceback = sys.exc_info()
			print "*** print_tb:"
			traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
			print "*** print_exception:"
			traceback.print_exception(exc_type, exc_value, exc_traceback,
			limit=2, file=sys.stdout)
			print "*** print_exc:"
			traceback.print_exc()
			print "*** format_exc, first and last line:"
			formatted_lines = traceback.format_exc().splitlines()
			print formatted_lines[0]
			print formatted_lines[-1]
			print "*** format_exception:"
			print repr(traceback.format_exception(exc_type, exc_value,
			exc_traceback))
			print "*** extract_tb:"
			print repr(traceback.extract_tb(exc_traceback))
			print "*** format_tb:"
			print repr(traceback.format_tb(exc_traceback))
			print "*** tb_lineno:", exc_traceback.tb_lineno

		finally:
			self.view.popupm.grab_release()

class TrajInfoView(tk.LabelFrame):
	'''Info'''
	def __init__(self,controller,parent=None):
		tk.LabelFrame.__init__(self,parent,text="Info:")
		self.controller=controller

		self.name = Line(self,"Name:")
		self.name.pack(side=tk.TOP, expand=1, fill=tk.X)
		self.ref = Line(self,"Reference:")
		self.ref.pack(side=tk.TOP, expand=1, fill=tk.X)
		self.frames = Line(self,"Frames:")
		self.frames.pack(side=tk.TOP, expand=1, fill=tk.X)
		self.atoms = Line(self,"Atoms:")
		self.atoms.pack(side=tk.TOP, expand=1, fill=tk.X)
		self.armsd = Line(self,"avg RMSD:")
		self.armsd.pack(side=tk.TOP, expand=1, fill=tk.X)

		#self.text.bind('<Button-2>',self.controller.onPopUp)
		#self.text.bind('<Button-3>',self.controller.onPopUp)


		self.figure = Plot.Figure(self)	

		self.makePopUp()
		if not self.name:
			print "hgskdgjkshdjkaJKSDjkah"
		self.update(None);

	def update(self,data):
		if data:
			self.name.data.set(data.name)
			self.ref.data.set(data.ref)
			self.frames.data.set(data.frames)
			self.atoms.data.set(data.atoms)
			self.armsd.data.set(data.armsd)
			self.figure.replot(range(15),range(15))
		else:
			self.name.data.set("None")
			self.ref.data.set("")
			self.frames.data.set("")
			self.atoms.data.set("")
			self.armsd.data.set("")

	def makePopUp(self):
		self.popupm = tk.Menu(self, tearoff=0) 
		self.popupm.add_command(label="Show")#,command=None)
				
