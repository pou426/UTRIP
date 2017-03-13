CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = mlwf.o mlwh.o mlwsh.o mlwsf.o mlwreset.o \
    getfro.o getfee.o getfrx.o getfbx.o getsaz.o getfeb.o \
    mlw_mh.o mlw_mh_c.o fx2c.o mlw_et.o mlw_gv.o mfdimto3.o \
    mlw_mh_s.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
