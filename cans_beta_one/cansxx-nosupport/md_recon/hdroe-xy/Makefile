CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o engine.o exp_et.o

$(EXE) : $(OBJS) $(CANS_LIB_DIR)/libcans3d.a
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
