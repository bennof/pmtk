'''

    ./python/ChiPlayer.py 

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





class ChiPlayerView(tk.LabelFrame):
	def __init__(self,controller,frame,root=None):
		tk.LabelFrame.__init__(self,parent,text="Player:")

            	button = tk.Button(self, text="Play", command=controller.onPlay).pack(side=tk.LEFT,padx=5)
            	button = tk.Button(self, text="Stop", command=controller.onStop).pack(side=tk.LEFT,padx=5)
            	button = tk.Button(self, text="Load", command=controller.onLoad).pack(side=tk.LEFT,padx=5)
		mode = tk.StringVar(self)
		option = tk.OptionMenu(self,mode, "Forward", "Backward", "Rock",command=controller.onModeChanged(mode))
		option.pack()
		if not root is None:
			root.bind("<p>", controller.onPlay)
			root.bind("<s>", controller.onStop)
			root.bind("<l>", controller.onLoad)

class ChiPlayer():
	def __init__(parent,frame):
		self.parent=parent
		self.view=PlayerView(self,frame)
		self.mode="Rock"
		self.view.mode.set(self.mode)
		self.trajs=[]
		self.trajsf=[]
		self.frames=0
		self.trigger=0
		self.cframe=0
	def onPlay(self):
            	if self.trigger:
                	CHIMERA.triggers.deleteHandler('new frame',self.trigger)
            	if self.frames:
                	self.trigger = CHIMERA.triggers.addHandler('new frame',self.playing,None)
	def onStop(self):
            if self.trigger:
                CHIMERA.triggers.deleteHandler('new frame',self.trigger)
                self.trigger=None
	def onLoad(self):
		obj = self.parent.getCur()
            	name   = obj.getName()
            	ref    = obj.references[name]
            	natoms = obj.getPoints();
            	frames = obj.getFrames()

		#get id of reference
            	for id in CHIMERA.openModels.list():
                	if id.name == ref:
                    		break;

            	if id.name == ref:
                	exframes=len(id.coordSets)
            	try: 
                	idx = self.trajs.index(id)
                	self.trajsf[idx]=frames
            	except Exception , e:
                	self.trajsf.append(frames)
                	self.trajs.append(id)
            
            	self.frames=frames
            	for j in self.trajf:
                	if j<self.frames:
                    		self.frames=j
            	at = id.atoms
		j=0
            	while j<exframes and j<frames:
                	cs=id.findCoordSet(j+1)
                	CHIMERA.fillCoordSet(cs,id.atoms,numpy.array(self.getCoords(self.parent.getFrame(j),natoms)))
                	j=j+1
            	while j<frames:
                	cs=id.newCoordSet(j+1,natoms)
                	CHIMERA.fillCoordSet(cs,id.atoms,numpy.array(self.getCoords(self.parent.getFrame(j),natoms)))
                	j=j+1

	def onModechange(self,mode):
		self.mode=mode
		pass
        def getCoords(self,cf,n):
            	j=0 
            	h=[]
            	while j < n:
                	x=cf[3*j]
                	y=cf[3*j+1]
                	z=cf[3*j+2]
                	j=j+1;
                	h.append(chi.Point(x,y,z))
            	return h
	
