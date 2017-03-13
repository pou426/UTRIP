CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans3d -lcansio

run : params.txt.0000
params.txt.0000 : $(EXE)
	mpirun -np 2 ./$(EXE)

clean :
	rm -f params.txt.* *.dac.* $(EXE) $(OBJS)

plot : idl.png
idl.png : params.txt.0000
	idl batch.pro
