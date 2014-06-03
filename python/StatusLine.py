'''

    ./python/StatusLine.py 

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

class StatusLine(tk.Label):
	'''Status line'''
	def __init__(self,parent):
		tk.Label.__init__(self,parent, bd=1, relief=tk.SUNKEN, anchor=tk.E)
		self.view=self
	def onUpdate(self,msg):
		self.config(text=msg)
		self.update_idletasks()
