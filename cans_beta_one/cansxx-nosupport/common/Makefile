CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = daupdate.o rv2vv.o bb2ex.o cfl_g.o cfl_e.o dx2dxi.o dx2ux.o \
vvtorv.o rvtovv.o ee2eh.o eh2pr.o ohmideal.o ohmet.o bbtocc.o \
pr2eh.o eh2ee.o adddfdx.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
