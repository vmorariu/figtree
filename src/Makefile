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
# If you do not have the ANN library on your system or wish to compile 
# the IFGT library without ANN support, either uncomment the line
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
# Set compiler and linker depending on whether the library is to be compiled
# as a shared library (*.so) or a static library (*.a).  The default is 
# a shared library.  However, if FIGTREE_LIB_TYPE=static then a static library 
# will be compiled.  This can be set either by uncommenting the line below, or
# by calling make with "FIGTREE_LIB_TYPE=static" as an argument:
#
#    make "FIGTREE_LIB_TYPE=static
#
#-----------------------------------------------------------------------------
#FIGTREE_LIB_TYPE=static

ifeq ($(FIGTREE_LIB_TYPE),static)
COMPILE=g++
MAKELIB=ar ruv
CFLAGS=-c -O3
LIBSUFFIX=.a
else
COMPILE=g++
MAKELIB=g++ -shared -o
CFLAGS=-c -O3 -fpic
LIBSUFFIX=.so
endif

CFLAGS+= -I$(INCDIR)

ifeq ($(FIGTREE_NO_ANN),true)
CFLAGS+= -DFIGTREE_NO_ANN
else
CFLAGS+= -I$(ANNINCDIR)
ifneq ($(FIGTREE_LIB_TYPE),static)
LDFLAGS= -L$(ANNLIBDIR) -lann_figtree_version
endif
endif

SOURCES=figtree.cpp KCenterClustering.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGETBASE=$(LIBDIR)/libfigtree
TARGET=$(TARGETBASE)$(LIBSUFFIX)

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(SOURCES) $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(MAKELIB) $@ $(OBJECTS) $(LDFLAGS)

.cpp.o:
	$(COMPILE) $(CFLAGS) $< -o $@

clean:
	rm -f ${TARGETBASE}.a ${TARGETBASE}.so ${OBJECTS}

mostlyclean:
	rm -f $(OBJECTS)
