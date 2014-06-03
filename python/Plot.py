'''

    ./python/Plot.py 

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

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from matplotlib.figure import Figure as matplotlibfigure

class Figure():
	def __init__(self,parent):
		self.f = matplotlibfigure(figsize=(5,4), dpi=100,facecolor='white')
		self.a = self.f.add_subplot(111)
		self.a.set_axis_bgcolor((1, 1, 1))
		self.a.plot(range(10),range(10))
		self.a.set_title('Wait ...')
		self.a.set_xlabel('X')
		self.a.set_ylabel('Y')
		self.canvas = FigureCanvasTkAgg(self.f, master=parent)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
	def replot(self,x,y,title="",xname="",yname=""):
		self.a.plot(x,y)
		self.a.set_title(title)
		self.a.set_xlabel(xname)
		self.a.set_ylabel(yname)
		self.canvas.show()
		pass

	
