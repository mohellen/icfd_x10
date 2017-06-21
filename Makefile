SRCDIR = src
BINDIR = bin
OUTDIR = out

ifndef MAINFILE
MAINFILE = $(SRCDIR)/Main.x10
endif

X10CPP = x10c++
X10CPPFLAGS =

all:
	mkdir -p $(OUTDIR)
	mkdir -p $(BINDIR)
	$(X10CPP) $(X10CPPFLAGS) $(MAINFILE) -d $(BINDIR) -o $(BINDIR)/Main

run:
	bin/Main

clean:
	rm -rf $(BINDIR) $(OUTDIR)
	rm -f *.trace

clear:
	rm -f $(OUTDIR)/* 
	rm -f *.trace
