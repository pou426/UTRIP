CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o bndcmat.o engine.o cndsor2.o ccfspt.o sorbr2.o gtnorm2.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans3d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)

plot : pldt.png
pldt.png : params.txt
	idl batch.pro
