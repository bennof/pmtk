'''

    ./python/CalcPanel.py 

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

class CalcPanel(tk.LabelFrame):
	'''Calculation lable '''
	def __init__(self,parent,controller):
		tk.LabelFrame.__init__(self,parent,text="Calculations:") 
		self.view=self
        	button = tk.Button(self, text="Align",   command=controller.onAlign).pack(side=tk.LEFT,padx=5)
        	button = tk.Button(self, text="Shift",   command=controller.onShift).pack(side=tk.LEFT,padx=5)
        	button = tk.Button(self, text="Average", command=controller.onAverage).pack(side=tk.LEFT,padx=5)
        	button = tk.Button(self, text="PCA",     command=controller.onPCA).pack(side=tk.LEFT,padx=5)
        	button = tk.Button(self, text="Make Traj.", command=controller.onTrajectory).pack(side=tk.LEFT,padx=5)
        	button = tk.Button(self, text="VRML",       command=controller.onVRML).pack(side=tk.LEFT,padx=5)
	
