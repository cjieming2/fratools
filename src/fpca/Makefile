DEBUG_OPTIONS= -g
NLIB=$(PWD)/INCLUDE/nicklib.a
IDIR=$(PWD)/INCLUDE
CFLAGS= -c -g -p -I$(IDIR) -Wimplicit-int

M1=fpca
M1O=fpca.o  eigsubs.o  eigx.o

$(M1): $(M1O)
	rm  -f  $(M1)
	gcc -static -I$(IDIR) $(DEBUG_OPTIONS) -o $(M1) $(M1O) ${NLIB} -lm -L${PWD} -llapack -lblas1 -lf2c

clean: 
	rm -f *.o 
	rm -f core

