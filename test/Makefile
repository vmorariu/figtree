#-----------------------------------------------------------------------------
# Some basic definitions:
#		BASEDIR		where include, src, lib, ... are
#		INCLIB		include directory
#		LIBLIB		library directory
#-----------------------------------------------------------------------------
BASEDIR	= ..
INCDIR	= $(BASEDIR)/include
LIBDIR	= $(BASEDIR)/lib
BINDIR      = $(BASEDIR)/bin

#-----------------------------------------------------------------------------
# Default locations for the Approximate Nearest Neighbor (ANN) library.
# For other locations, modify the paths in this makefile, or set them on the
# command line.  For example, if the ANN library resides in your_ann_path,
# you can either modify the variables below, or call 'make' as below:
#
#    make "ANNLIBDIR=your_ann_path/lib" "ANNINCDIR=your_ann_path/include"
#
# If the FIGTree library was compiled without ANN support, either uncomment the line
# that defines "FIGTREE_NO_ANN=true" or set it on the command line:
#
#    make "FIGTREE_NO_ANN=true"
#
#-----------------------------------------------------------------------------
ANNLIBDIR   = $(BASEDIR)/lib
ANNINCDIR   = $(BASEDIR)/external/ann_1.1.1/include

# uncomment this line to compile the library without ANN support
#FIGTREE_NO_ANN=true

#-----------------------------------------------------------------------------
# Set compiler and linker, as well as their flags, sources, and targets.
#-----------------------------------------------------------------------------
CC=g++
CFLAGS=-c -O3 -I$(INCDIR)
LDFLAGS=-L$(LIBDIR) -lfigtree

# if the FIGTree library was compiled with ANN support, then link with ANN library
# (this is really necessary only if ANN and FIGTree are compiled as static libs,
# but doesn't hurt in the other cases)
ifneq ($(FIGTREE_NO_ANN),true)
LDFLAGS+= -L$(ANNLIBDIR) -lann_figtree_version
endif

SOURCES=test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=$(BINDIR)/test


#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(SOURCES) $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $< $(CFLAGS) -o $@

clean:
	rm -f ${TARGET} ${OBJECTS}

mostlyclean:
	rm -f $(OBJECTS)

