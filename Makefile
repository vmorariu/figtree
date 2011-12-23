all clean mostlyclean:
	cd src; $(MAKE) $@
	cd samples; $(MAKE) $@
	cd test; $(MAKE) $@

