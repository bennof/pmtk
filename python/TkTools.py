'''

    ./python/TkTools.py 

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
import tkFileDialog as tkFile
import tkMessageBox as tkMessage


def openDialog(where,filetypes=[('all files', '.*')],multiple=True,title='Open Trajectory',message='Open Trajectory (PDB, XTC)',):
	return  tkFile.askopenfilename(title=title,message=message,initialdir=where,multiple=multiple)
def safeDialog(where):
	return tkFile.askopenfilename(parent=self,initialdir=where)
def warning(msg):
	tkMessage.showwarning("Warning","Warning: "+msg)

class Line(tk.Frame):
	def __init__(self,parent,name):
		tk.Frame.__init__(self,parent)
		self.data = tk.StringVar()
		tk.Label(self, text=name).pack(side=tk.LEFT);
		tk.Label(self,textvariable=self.data).pack(side=tk.LEFT);
		pass;
	
class Entry(tk.Frame):
	def __init__(self,parent,name):
		tk.Frame.__init__(self,parent)
		self.data = tk.StringVar()
		tk.Label(self, text=name).pack(side=tk.LEFT);
		tk.Entry(self,textvariable=self.data).pack(side=tk.LEFT);
		pass;

class DialogRef(tk.Toplevel):
	def __init__(self,parent,controller,list ={}):
		tk.Toplevel.__init__(self,parent)
		self.parent=parent;
		self.title("Select Reference")
		tk.Label(self, text="Models:").pack(side=tk.TOP)
		self.List = tk.Listbox(self,selectmode=tk.SINGLE)
		self.List.pack(fill=tk.BOTH, expand=1)

            	for i in list:
                	self.List.insert(tk.END, i.name)

            	frame = tk.Frame(self)
            	button = tk.Button(frame, text="select", command=self.onSelect)
            	button.pack(side=tk.LEFT)
            	button = tk.Button(frame, text="open", command=self.onOpen)
            	button.pack(side=tk.RIGHT)
            	button = tk.Button(frame, text="close", command=self.onExit)
            	button.pack(side=tk.RIGHT)
            	frame.pack(side=tk.BOTTOM)  

		self.List.bind('<Double-Button-1>', self.onSelect)

		self.grab.set()
		self.wait_window(self)
	
	def onOpen(self):
		pass
	def onSelect(self):
		pass
	def onExit(self):
		self.parent.focus_set()
		self.destroy()
