TARGET          = xdrfile
FOLDER          = xdrfile-1.1b
BUILDDIR        = $(FOLDER)
TARFILE         = xdrfile-1.1.tar.gz
URL             = ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz 
CONFIGURE_FLAGS = 

all: $(TARGET).maked
	-echo "xdr/xdrfile-1.1b/src/.libs/libxdrfile.a"
      
install:

getstaticlibpath:
	@echo "$(shell pwd)/xdrfile-1.1b/src/.libs/libxdrfile.a"

	 
$(TARGET).installed : $(TARGET).maked  
	$(MAKE) -C $(BUILDDIR) install

$(TARGET).maked : $(TARGET).configured
	$(MAKE) -C $(BUILDDIR)
	-echo "$(BUILDDIR): make " > $@
 
$(TARGET).configured : $(FOLDER)
	echo $(LDFLAGS)
	-mkdir -p $(BUILDDIR)
	cd $(BUILDDIR) && ../$(FOLDER)/configure $(CONFIGURE_FLAGS)
	-echo "$(BUILDDIR): configure $(CONFIGURE_FLAGS)" > $@

$(FOLDER) : $(TARFILE)
	$(TAR) zxf $(TARFILE) || $(TAR) jxf  $(TARFILE) || $(TAR) Jxf $(TARFILE)

$(TARFILE) :
	$(DOWNLOADER) $(DOWNLOADER_ARGS) $@ $(URL) 

clean:
	-rm -rf $(FOLDER) $(TARFILE) $(BUILDDIR) $(TARGET).configured $(TARGET).maked  $(TARGET).installed
