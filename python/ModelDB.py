'''

    ./python/ModelDB.py 

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



import os
class Model:
	def __init__(self,name):
		self.name=name
		self.ref=None
		self.frames=None
		self.atoms=None
		self.armsd=None


class ModelDB:
	def __init__(self):
		self.data={}
	def new(self,name,fn=None,model=None):
		if fn:
			pass
		traj = Model(name)
			#traj.read(fn)
		if name in self.data:
			self.data[name+"_new"]=traj
		else:
			self.data[name]=traj
		return traj;
	def delete(self,name):
		if name in self.data:
			del self.data[name]
		return None
	def get(self,name):
		if name in self.data: 
			return self.data[name];
		else:
			return None

	def getNames(self):
		return self.data.keys()
