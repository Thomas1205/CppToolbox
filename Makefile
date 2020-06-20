include ../Makefile.common

all: $(BIN)/testcommon.debug $(BIN)/testcommon.opt  $(LIB)/commonlib.debug $(LIB)/commonlib.opt $(DEBUGDIR)/vector.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o 

$(BIN)/testcommon.debug: test.cc matrix.hh trimatrix.hh storage1D.hh storage2D.hh tristorage2D.hh storage3D.hh vector.hh application.hh makros.hh routines.hh sorting.hh grayimage.hh tensor.hh colorimage.hh fileio.hh matrix_inversion.hh timing.hh outer_product.hh hash_map.hh unsorted_map.hh sorted_map.hh unsorted_set.hh sorted_set.hh $(DEBUGDIR)/timing.o $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o $(DEBUGDIR)/combinatoric.o $(DEBUGDIR)/rational.o $(DEBUGDIR)/makros.o
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) test.cc $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(DEBUGDIR)/combinatoric.o $(DEBUGDIR)/rational.o $(DEBUGDIR)/makros.o -o $@

$(BIN)/testcommon.opt: test.cc matrix.hh trimatrix.hh storage1D.hh storage2D.hh tristorage2D.hh storage3D.hh vector.hh application.hh makros.hh routines.hh grayimage.hh tensor.hh colorimage.hh fileio.hh matrix_inversion.hh timing.hh $(OPTDIR)/timing.o $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o $(OPTDIR)/combinatoric.o $(OPTDIR)/rational.o $(OPTDIR)/makros.o
	$(LINKER) $(OPTFLAGS) $(INCLUDE) test.cc $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(DEBUGDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/combinatoric.o $(OPTDIR)/rational.o $(OPTDIR)/makros.o -o $@

$(LIB)/commonlib.debug: $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(DEBUGDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o $(DEBUGDIR)/makros.o $(DEBUGDIR)/combinatoric.o
	ar rs $@ $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(DEBUGDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o  $(DEBUGDIR)/makros.o  $(DEBUGDIR)/combinatoric.o
	ranlib $@

$(LIB)/commonlib.opt: $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o 
	ar rs $@ $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o
	ranlib $@

$(BIN)/testasm.opt.L64: testasm.cc $(DEBUGDIR)/timing.o
	$(LINKER) $(DEBUGFLAGS) $(INCLUDE) testasm.cc $(DEBUGDIR)/timing.o -o $@

include ../Makefile.finish

