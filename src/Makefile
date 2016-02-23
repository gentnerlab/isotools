
CC                     = g++
CFLAGS = -Wall -O2

BINDIR = .
INSTALL = install
INSTALLBIN = -s -m 755
INSTALLMAN = -m 644
SUBDIRS = isoi isorat wave2f
#MYMAKEFLAGS = 'BINDIR=$(BINDIR)' 'MANDIR=$(MANDIR)' 'MANSECT=$(MANSECT)'\
# 'INSTALL=$(INSTALL)' 'INSTALLBIN=$(INSTALLBIN)' 'INSTALLMAN=$(INSTALLMAN)'\
# 'CC=$(CC)' 'SWITCHES=$(SWITCHES)' 'O=$(O)' 'LIBRARYOBJS=$(LIBRARYOBJS)' \
# 'OPTIMIZE=$(OPTIMIZE)' 'MACHINE=$(MACHINE)'

all:
	@for i in $(SUBDIRS); do \
	echo "make all in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) all); done

install:
	@for i in $(SUBDIRS); do \
	echo "Installing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS)  install); done

clean:
	@for i in $(SUBDIRS); do \
	echo "Clearing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) clean); done
