CANS_LIB_DIR = ../..

OBJS = main.o model.o bnd.o cipbnd.o 

params.txt : a.out
	./a.out

a.out : $(OBJS)  $(CANS_LIB_DIR)/libcans1d.a
	$(FC) -o a.out $(OBJS) \
	-L$(CANS_LIB_DIR) -lcans1d -lcansnc

clean :
	rm -f params.txt *.dac a.out $(OBJS) 
