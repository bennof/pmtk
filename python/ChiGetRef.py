'''

    ./python/ChiGetRef.py 

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



import TKinter as yk

class  ChiGetRefView(tk.Toplevel):
	''' Creates a window to select a reference in chimera (MVC - Model)'''
	def __init__(self,controller,parent = None,list = {}):
		''' init '''
		tk.Toplevel.__init__(self,parent) 
		self.parent=parent
            	self.title("Select Reference")
            	self.protocol("WM_DELETE_WINDOW",controller.onExit)
            	tk.Label(self, text="Models:").pack(side=tk.TOP)

                self.List = tk.Listbox(self,selectmode=tk.SINGLE)
                self.List.pack(fill=tk.BOTH, expand=1)
		
            	for i in list:
                	self.List.insert(tk.END, i.name)

            	frame = tk.Frame(self)
            	button = tk.Button(frame, text="select", command=controller.onSelect)
            	button.pack(side=tk.LEFT)
            	button = tk.Button(frame, text="open", command=controller.onOpen)
            	button.pack(side=tk.RIGHT)
            	frame.pack(side=tk.BOTTOM)  
            	self.chiList.bind('<Double-Button-1>', controller.onSelect)
	def getSelect():
		''' get current selection '''
		return self.List.curselection() 
	def update(List):
		for i in List: 
			delete(0, END) 
			self.List.insert(tk.END, i.name)

class ChiGetRef():
	''' controller for selcting a reference in chimera '''
	def __init__(self,parent=None):
		''' init '''
	 	self.controller=parent
		self.view=ChiGetRefView(self,parent.view)
		pass
	def getRef(self):
	     	''' gets data from Cimera opens view'''   
		if not self.view is None:
			self.view = ChiGetRefView(self,CHIMERA.openModels.list())
	def onExit(self):
		''' exit event '''
		self.view.destroy()
		self.view=None
	def onOpen(self):
		''' open event '''
	   	cwd = None
		self.controlle.sendStatus("Open Reference ...")
		if not self.controller is None:
			cwd = getCWD();
		name = openDialog(cwd)
		if name:
			self.controller.sendStatus("Open Reference: "+name+" ...")
			id = CHIMERA.openModels.open(name)
			id=len(CHIMERA.openModels.list())-1
			obj = CHIMERA.openModels.list()[id]
			self.controller.sendStatus("Checking Reference: "+name+" ...")
			if len(obj.atoms)!=self.controller.getCurrent().getAtoms():


				warning("Refrence differs in Structure")
				self.controller.sendStatus("Checking Reference: "+name+" Failed")
				self.controller.setCurRef(None);
			else:
				self.controller.sendStatus("Checking Reference: "+name+" ok")
				self.controller.setCurRef(i);
		else:
			self.controller.sendStatus("Open Reference: None")
			self.controller.setCurRef(None);
	def onClose(self):
		''' close model in chimera event '''
		id = self.view.getSelect()	
		r = CHIMERA.openModels.list()[int(id[0])]
		self.controller.sendStatus("Close Reference: "+r.name)
		CHIMERA.openModels.close(r)
		self.controller.sendStatus("Update References")
		self.view.update(CHIMERA.openModels.list())
	def onSelect(self,):
		id = self.view.getSelect()	
		r = CHIMERA.openModels.list()[int(id[0])]
		self.controller.sendStatus("Selecting: "+r.name )
		self.controller.setCurRef(r);
		onExit()

	
