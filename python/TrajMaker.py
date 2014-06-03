'''

    ./python/TrajMaker.py 

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

class  TrajMakerView(tk.Toplevel):
	def __init__(self,controller,parent = None,list = {}):
		''' init '''
		tk.Toplevel.__init__(self,parent) 
		self.parent=parent
		self.title("Make Trajectory")
            	self.protocol("WM_DELETE_WINDOW",controller.onExit)

		self.namev   = tk.StringVar()
		self.namev.set("NewTrajectory")
		self.namel   = tk.Label(self, text="Name:")
		self.namee   = tk.Entry(self,textvariable=self.namev)
            	#self.nameb   = tk.Button(self, text="new", command=controller.onNew)
		self.namel.grid(row=0, column=0,sticky=tk.W)
		self.namee.grid(row=0, column=1,sticky=tk.W)
		#self.nameb.grid(row=0, column=2,sticky=tk.W)
		
		self.refv   = tk.StringVar()
		self.refl   = tk.Label(self, text="Reference:")
		self.refe   = tk.Entry(self,textvariable=self.refv)
            	self.refb   = tk.Button(self, text="select", command=controller.onSelectRef)
                self.refl.grid(row=1, column=0,sticky=tk.W)
                self.refe.grid(row=1, column=1,sticky=tk.W)
                self.refb.grid(row=1, column=2,sticky=tk.W)

		self.vecv    = tk.StringVar()
		self.vecl    = tk.Label(self, text="Vector:")
		self.vece    = tk.Entry(self,textvariable=self.vecv)
            	self.vecb    = tk.Button(self, text="select", command=controller.onSelectVec)
		self.vecfv   = tk.StringVar()
		self.vecfl   = tk.Label(self, text="Frame:")
		self.vecfe   = tk.Entry(self,textvariable=self.vecfv,width=6)
                self.vecl.grid(row=2, column=0,sticky=tk.W)
                self.vece.grid(row=2, column=1,sticky=tk.W)
                self.vecb.grid(row=2, column=2,sticky=tk.W)
                self.vecfl.grid(row=2, column=3,sticky=tk.W)
                self.vecfe.grid(row=2, column=4,sticky=tk.W)

		self.lengthv = tk.StringVar()
		self.lengthv.set(10.0);
		self.lengthl = tk.Label(self, text="Length in (A):")
		self.lengthe = tk.Entry(self,textvariable=self.lengthv,width=10);
                self.lengthl.grid(row=3, column=0,sticky=tk.W)
                self.lengthe.grid(row=3, column=1,sticky=tk.W)

		self.framesv = tk.StringVar()
		self.framesv.set(10)
		self.framesl = tk.Label(self, text="Frames:")
		self.framese = tk.Entry(self,textvariable=self.framesv,width=10);
                self.framesl.grid(row=3, column=2,sticky=tk.W)
                self.framese.grid(row=3, column=3,sticky=tk.W)

		ff = tk.Frame(self)
		ff.grid(row=4, column=0,sticky=tk.W,rowspan=6)
		make    = tk.Button(ff, text="Make", command=controller.onMake)		
		#preview = tk.Button(self, text="Make", command=controller.onMake)		
		close   = tk.Button(ff, text="Close", command=controller.onExit)		
		make.pack(side=tk.RIGHT)
		close.pack(side=tk.RIGHT)
		#self.minsize(300,300)

class TrajMaker():
	def __init__(self,controller):
		self.controller=controller
		self.view = TrajMakerView(self,controller.view)
	def onExit(self):
		self.view.destroy()
		self.view=None
	def onSelectVec(self):

		pass
	def onSelectRef(self):
		DialogRef(self.view,self)	
		pass
	def onMake(self):
		pass
	def onNew(self):
		warning("Method removed")	
		pass
		

