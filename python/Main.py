'''

    ./python/Main.py 

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



from ModelDB import *
from CalcPanel import *
from StatusLine import * 
from TrajList import *
from TrajInfo import *
from TrajMaker import *
from TkTools import *


class Main():
	def __init__(self,name="Test work"):
		#self init
		self.cwd=""
		
		#Chimera connection

		#init data
		self.data=ModelDB()
		self.cur=None

		#init view
		self.view=MainView(self,name)


		self.view.select.activate()
		self.updateAll();
		self.view.status.onUpdate("Done")

	def run(self):
		self.view.mainloop()
	def getCWD(self):
		return self.cwd
	def sendStatus(self,msg):
		self.view.status.onUpdate(msg)
	def updateAll(self):
		self.view.info.onUpdate(self.cur)
		self.view.select.onUpdate(self.data.getNames(),self.cur)
	def onExit(self):
		self.view.destroy()
	def onNew(self):
		fname = saveDialog(self.cwd)	
		self.cur=self.data.new(fname)		
		self.view.select.update(self.getTrajList())
		pass
	def onChangeCurrent(self):
		idx = self.view.select.List.curselection() 
		print idx;
	def onOpen(self,fnames=None):
		self.sendStatus("Open ...")
	  	if not fnames:
			fnames = openDialog(self.cwd)	
		if not fnames:
			self.sendStatus("Open None")
			return
		for fname in fnames:
			if fname:
				self.sendStatus("Open "+fname+" ...")
				self.cur=self.data.new(os.path.basename(fname),fname)		
				self.updateAll();
		self.sendStatus("Done")
		pass
	def onSave(self,fname=None):
	        if not self.cur:
			warning("No date selected to save")
			self.view.status.onUpdate("Save None")
		self.view.status.onUpdate("Save ...")
	  	if not fname:
			fname = saveDialog(self.cwd)	
		if not fname:
			self.view.status.onUpdate("Save None")
		self.view.status.onUpdate("Save "+fname+" ...")
		self.view.select.update(self.getTrajList())
		pass
	def onAppend(self,name=None,fielen=None):
		self.sendStatus("Append ...")
	  	if not fname:
			fname = openDialog(self.cwd)	
		if not fname:
			self.view.status.onUpdate("Appand None")
			return
		self.view.select.onUpdate(self.data.getNames())
		pass
	def onDelete(self,name=None):
		if not name:
			name=self.cur.name;
			self.data.delete(name)
			self.cur=None
		else:
			self.model.delete(name)
		self.updateAll();
	def onAlign(self):
		pass
	def onShift(self):
		pass
	def onAverage(self):
		pass
	def onPCA(self):
		pass
	def onTrajectory(self):
		TrajMaker(self)
		pass
	def onVRML(self):
		pass


class MainView(tk.Tk):
	def __init__(self,controller,name=None,parent=None,player=False):
                tk.Tk.__init__(self,parent)
		print __file__
		self.parent=parent
		self.controller=controller

		#tk.Toplevel.__init__(self, parent, width=600,height=500) 

            	self.grid_propagate(0)
            	self.title(name)
            	self.protocol("WM_DELETE_WINDOW", self.controller.onExit)
		self.minsize(300,300)
		self.geometry("600x500")

		self.makeMenu()
		self.makeFrame(player)

	def makeMenu(self):
		""" creates main menu """
            	menubar = tk.Menu(self,tearoff=0)
		#File Menu
		file = tk.Menu(menubar, tearoff=0) 
                file.add_command(    label="Open Trajectory",     command=self.controller.onOpen, accelerator="Ctrl+o")
            	file.add_command(    label="Append Trajectory",   command=self.controller.onAppend, accelerator="Ctrl+a")
            	file.add_command(    label="Save Trajectory",     command=self.controller.onSave, accelerator="Ctrl+s")
            	file.add_command(    label="Delete Trajectory",   command=self.controller.onDelete, accelerator="Ctrl+d")
                file.add_separator()
                file.add_command(    label="Quit",                command=self.controller.onExit,           accelerator="Ctrl+q")
                menubar.add_cascade(label="File", menu=file)

            	#Calculate Menu
            	calc = tk.Menu(menubar, tearoff=0)
            	calc.add_command(label="Align",     command=self.controller.onAlign,      accelerator="Alt+A")
            	calc.add_command(label="Shift",     command=self.controller.onShift,      accelerator="Alt+S")
            	calc.add_command(label="Average",   command=self.controller.onAverage,    accelerator="Alt+M")
            	calc.add_command(label="PCA",       command=self.controller.onPCA,        accelerator="Alt+P")
            	calc.add_command(label="Trajectory",command=self.controller.onTrajectory, accelerator="Alt+T")
            	calc.add_command(label="VRML",      command=self.controller.onVRML,       accelerator="Alt+V")
            	menubar.add_cascade(label="Calculate", menu=calc)
		#Add menu
            	self.config(menu=menubar)

		pass
	def makeFrame(self,player):
		""" creates the frame """

		tk.Label(self, text="To load trajectories into Chimera please select a reference.").pack(side=tk.TOP)   

            	self.status = StatusLine(self)
            	self.calc = CalcPanel(self,self.controller)
		hframe=tk.Frame(self);
		self.select = TrajList(hframe,self.controller)
		self.info = TrajInfo(hframe)

		if player is True: 
		        pc = PlayerController()
            		self.player=PlayerView(self,pc,root=self)
			self.player.pack(side=tk.BOTTOM,fill=tk.X)


		self.status.view.pack(fill=tk.X,side=tk.BOTTOM)
		self.calc.view.pack(side=tk.BOTTOM,fill=tk.X)
		self.select.view.pack(side=tk.LEFT,fill=tk.BOTH,expand=1)
		self.info.view.pack(side=tk.RIGHT,fill=tk.BOTH,expand=1)
		hframe.pack(side=tk.TOP,fill=tk.BOTH,expand=1)

	def status(self,msg):
		self.statusline.set(msg);
		
	


