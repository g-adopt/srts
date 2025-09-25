#-- Gfortran

FCOMPILE  = gfortran -c  -ffixed-line-length-none         -fno-automatic -O3

MAKE          = make

.f.a:
	${FCOMPILE}  $<
	ar rv $@ $*.o
	rm -f $*.o
.F.a:
	${FCOMPILE}  $<
	ar rv $@ $*.o
	rm -f $*.o

LIB = libSphFunc
LIBDIR = ..

$(LIBDIR)/$(LIB).a: \
	$(LIBDIR)/$(LIB).a(legndr.o) \
	$(LIBDIR)/$(LIB).a(normylm.o) \
	$(LIBDIR)/$(LIB).a(ylm.o) \
	$(LIBDIR)/$(LIB).a(ylm02xlm.o)

	ranlib $(LIBDIR)/$(LIB).a

# Clean target
.PHONY: clean
clean:
	rm -f *.o $(LIBDIR)/$(LIB).a
