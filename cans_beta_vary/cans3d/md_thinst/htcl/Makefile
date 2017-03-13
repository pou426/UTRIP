CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o pertub.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
                -L$(CANS_LIB_DIR) -lcans3d -lcansio

run : params.txt
params.txt : a.out
	./a.out

clean:
	rm -f params.txt *.dac $(EXE) $(OBJS)

plot : pldt.png
pldt.png : params.txt
	idl batch.pro
