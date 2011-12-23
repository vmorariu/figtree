ifeq ($(FIGTREE_LIB_TYPE),static)
LIBSUFFIX=.a
MAKEOPT=linux-g++
else
LIBSUFFIX=.so
MAKEOPT=linux-g++-sl
endif

all:
	cd ann_1.1.1 ; $(MAKE) $(MAKEOPT)
	mv ann_1.1.1/lib/libann_figtree_version$(LIBSUFFIX) ../lib/libann_figtree_version$(LIBSUFFIX)

clean:
	cd ann_1.1.1 ; $(MAKE) realclean
	rm -f ../lib/libann_figtree_version.so
	rm -f ../lib/libann_figtree_version.a

mostlyclean:
	cd ann_1.1.1 ; $(MAKE) realclean
	rm -f ../lib/libann_figtree_version.so
	rm -f ../lib/libann_figtree_version.a

