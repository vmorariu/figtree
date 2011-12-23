all clean mostlyclean:
	cd external; $(MAKE) $@
	cd src; $(MAKE) $@
	cd samples; $(MAKE) $@
	cd test; $(MAKE) $@
