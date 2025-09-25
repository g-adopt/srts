FCOMPILE = gfortran -c -ffixed-line-length-none -fno-automatic -O3


MAKE = make
LIBDIR = ..

.f.a:
	${FCOMPILE}  $<
	ar rv $@ $*.o
	rm -f $*.o

$(LIBDIR)/libS20.a: \
	$(LIBDIR)/libS20.a(addata.o) \
	$(LIBDIR)/libS20.a(convpack.o) \
	$(LIBDIR)/libS20.a(dcopy.o) \
	$(LIBDIR)/libS20.a(scopy.o) \
	$(LIBDIR)/libS20.a(dot.o) \
	$(LIBDIR)/libS20.a(getfdp.o) \
	$(LIBDIR)/libS20.a(integ.o) \
	$(LIBDIR)/libS20.a(istlen.o) \
	$(LIBDIR)/libS20.a(inittrf.o) \
	$(LIBDIR)/libS20.a(legndr.o) \
	$(LIBDIR)/libS20.a(modtrf.o) \
	$(LIBDIR)/libS20.a(rmod.o) \
	$(LIBDIR)/libS20.a(saxpy.o) \
	$(LIBDIR)/libS20.a(sdot.o) \
	$(LIBDIR)/libS20.a(splh.o) \
	$(LIBDIR)/libS20.a(strgrep.o) \
	$(LIBDIR)/libS20.a(rsphhead.o) \
	$(LIBDIR)/libS20.a(rsple.o) \
	$(LIBDIR)/libS20.a(rspln.o) \
	$(LIBDIR)/libS20.a(splhsetup.o) \
	$(LIBDIR)/libS20.a(wsphhead.o) \
	$(LIBDIR)/libS20.a(wspthead.o) \
	$(LIBDIR)/libS20.a(wint2ch.o) \
	$(LIBDIR)/libS20.a(ylm.o)

	ranlib $(LIBDIR)/libS20.a

# Clean target
.PHONY: clean
clean:
	rm -f *.o $(LIBDIR)/libS20.a
