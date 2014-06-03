pmtk_PYTHON = 
   CalcPanel.py \
   ChiGetRef.py \
   ChiPlayer.py \
   Main.py \
   ModelDB.py \
   Plot.py \
   StatusLine.py \
   TkTools.py \
   TrajInfo.py \
   TrajList.py \
   TrajMaker.py \
   __init__.py 

pmtkdir = $(pythondir)/python

bin_SCRIPTS = pmtk
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = pmtk.in

do_substitution = sed -e 's,[@]pythondir[@],$(pythondir),g' \
	              -e 's,[@]PACKAGE[@],$(PACKAGE),g' \
		      -e 's,[@]VERSION[@],$(VERSION),g' 



pmtk : pmtk.in
	$(do_substitution) < $(srcdir)/pmtk.in > pmtk
	chmod +x pmtk




hellodir = $(pythondir)/hello


