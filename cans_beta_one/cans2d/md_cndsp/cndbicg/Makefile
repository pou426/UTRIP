CANS_LIB_DIR = ../../..

OBJS = model.o bnd.o bndcmat.o

a.out : main.o $(OBJS)
	$(FC) -o a.out main.o $(OBJS) \
	-L$(CANS_LIB_DIR) -lcans2d -lcansio

run : params.txt
params.txt : a.out
	./a.out

clean :
	rm -f params.txt *.dac a.out *.o

plot : pldt.png
pldt.png : params.txt
	idl batch.pro
