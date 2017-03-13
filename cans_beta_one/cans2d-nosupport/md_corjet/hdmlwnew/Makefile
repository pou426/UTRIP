CANS_LIB_DIR = ../../..
EXE = a.out

OBJS = main.o model.o bnd.o pertub.o etaanom.o \
  mlw_mh.o avlap_6.o rv2vv.o ee2pr.o diffcen.o bb2ez.o \
  getfro.o getfee.o getfrx.o getfry.o getfbx.o getfby.o getsaz.o \
  mlwh.o mlwf.o mlwreset.o daupdate.o mlw_g.o mlw_e.o mlwsf.o mlwsh.o cfl_e.o cfl_g.o \
  dx2dxi.o dx2ux.o engine.o

$(EXE) : $(OBJS)
	$(FC) -o $(EXE) $(OBJS) \
		-L$(CANS_LIB_DIR) -lcans2d -lcansio

run : params.txt
params.txt : $(EXE)
	./$(EXE)

clean :
	rm -f params.txt *.dac $(EXE) $(OBJS)

plot : pldt.png
pldt.png : params.txt
	idl batch.pro
