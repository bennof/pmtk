'''

    ./python/TrajList.py 

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

class TrajList(tk.LabelFrame):
	'''Trajectory list'''
	def __init__(self,parent,controller):
		tk.LabelFrame.__init__(self,parent,text="Trajectories:")
		self.List = tk.Listbox(self,bd=0)
		self.List.pack(side=tk.LEFT,fill=tk.BOTH, expand=1)
		frame=tk.Frame(self);
        	button = tk.Button(frame, text="+",   command=controller.onOpen).pack(side=tk.LEFT,padx=1)
        	button = tk.Button(frame, text="-",   command=controller.onDelete).pack(side=tk.LEFT,padx=1)
        	button = tk.Button(frame, text="s",   command=controller.onSave).pack(side=tk.LEFT,padx=1)
        	button = tk.Button(frame, text="++",  command=controller.onAppend).pack(side=tk.LEFT,padx=1)
		frame.pack(side=tk.BOTTOM,)

		self.controller=controller
		self.view=self
		self.config(width=200)
	def activate(self):
		self.List.bind('<ListboxSelect>',self.controller.onChangeCurrent())
	def onUpdate(self,List,cur):
		self.List.delete(0, tk.END)
		i=0;
		id=-1
		for item in List:

			if item == cur:
				id=i
			self.List.insert(tk.END,item)
			i=i+1
		if cur and id>=0:
			self.List.activate(id)
