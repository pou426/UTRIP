CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o bndsc.o cipbnd.o cip_mt_bg.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans1d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)

plot : pldt.png
pldt.png : params.txt
	idl batch.pro
