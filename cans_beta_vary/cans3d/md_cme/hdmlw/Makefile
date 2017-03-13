CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o mlw_m_sg2.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans3d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)

plot : idl.png
idl.png : params.txt
	idl batch.pro
