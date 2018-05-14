
SUBDIRS = log hnk oset sphericalwavelet tracking wavelet wavetree

default : all

all :
	$(foreach subdir,$(SUBDIRS),$(MAKE) -C $(subdir) all;) 

clean :
	$(foreach subdir,$(SUBDIRS),$(MAKE) -C $(subdir) clean;)

BASESOURCES=LICENSE Makefile README.md
DIR = TDTbase
dist :
	mkdir -p $(DIR)
	$(foreach subdir,$(SUBDIRS),$(MAKE) -C $(subdir) dist && tar -C $(DIR) -xzf $(subdir)/$(subdir).tar.gz;)
	$(foreach src,$(BASESOURCES),cp $(src) $(DIR);)
	tar -czf $(DIR).tar.gz $(DIR)/*
	rm -rf $(DIR)

