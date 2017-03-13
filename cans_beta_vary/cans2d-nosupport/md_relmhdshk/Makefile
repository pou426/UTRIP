CANS_LIB_DIR = ../..
OBJS = main.o model.o bnd.o

params.txt : a.out
	./a.out

a.out : $(OBJS) $(CANS_LIB_DIR)/libcans2d.a
	$(FC) $(FFLAGS) -o a.out $(OBJS) \
	-L$(CANS_LIB_DIR) -lcans2d -lcansnc

clean:
	rm -f params.txt *.dac a.out $(OBJS)
