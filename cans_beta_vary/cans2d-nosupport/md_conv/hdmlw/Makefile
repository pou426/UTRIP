CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o pertub.o mlw_h_tfix.o visten.o cf_tfxh.o cf_tfxf.o cfl_h_kv.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans2d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)
