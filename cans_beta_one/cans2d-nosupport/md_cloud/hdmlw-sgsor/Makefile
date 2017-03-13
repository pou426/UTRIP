CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o bndcmat.o bmperx.o bmpery.o selfgsor.o selfgbicg.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans2d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)
