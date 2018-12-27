CC=/usr/bin/mpicc
CFLAGS= -Wall -O3 -c
LD_FLAGS=-lgsl -lgslcblas -lm
GLIB_LOC = pkg-config --cflags --libs glib-2.0
GLIB_FLAGS = pkg-config --libs glib-2.0
TEX=/usr/bin/pdflatex
#Note: To disable the use of ATLAS for BLAS operations, comment the line below
#LD_FLAGS=-lgsl -lcblas -latlas -lm

WORK=${PWD}
CC_COMPILE=$(CC) `$(GLIB_LOC)` $(CFLAGS) $< 
TAR=/bin/tar
NAME=fermion_rand_pdrive
BIN_NAMES=isingrand_parallel
DIST=Makefile *.c *.h scripts README.tex
DOC=README.tex

#Dependencies
OBJECTS=initconds.o responses.o gsl_determinant_complex.o cantor_pair.o

.DEFAULT_GOAL := Release

Release: $(BIN_NAMES)

indent:
	indent *.c *.h
	${RM} *~

%.o:%.c %.h params.h
	$(CC_COMPILE)

isingrand_parallel: isingrand_parallel.o integrator_isingrand.o $(OBJECTS) params.h
	$(CC) isingrand_parallel.o integrator_isingrand.o $(OBJECTS) `$(GLIB_FLAGS)` $(LD_FLAGS) -o $@

dist:
	$(TAR) czvf $(NAME).tar.gz $(DIST)
	
doc:	
	$(TEX) $(DOC)
	$(TEX) $(DOC)

allclean: clean resclean

clean:
	${RM} $(BIN_NAMES) result/*.dat* *~ *.o *out *.dat* *.tgz *.aux *.log *.pdf *.synctex.gz

resclean:
	${RM} *result/*
	${RM} -r *result_*
