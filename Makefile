include Makefile.common

all: $(LIB)/commonlib.debug $(LIB)/commonlib.opt $(DEBUGDIR)/vector.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o 

$(LIB)/commonlib.debug: $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(DEBUGDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o $(DEBUGDIR)/makros.o $(DEBUGDIR)/combinatoric.o
	ar rs $@ $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(DEBUGDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o  $(DEBUGDIR)/makros.o  $(DEBUGDIR)/combinatoric.o
	ranlib $@

$(LIB)/commonlib.opt: $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o 
	ar rs $@ $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o
	ranlib $@

include Makefile.finish

