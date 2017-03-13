CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans3d -lcansio

run : out.txt
out.txt : $(EXE)
	./$(EXE)

clean :
	rm -f out.txt params.txt *.dac $(EXE) $(OBJS)

plot : pldt.png
pldt.png : out.txt
	idl batch.pro
