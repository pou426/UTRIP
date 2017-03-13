CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = gtnorm.o cndprm.o \
	bicgstab1.o bicgstab2.o iludcmp.o psolv.o \
	cndbicg.o ccfspt.o \
	cndbicg_c.o ccfspt_c.o \
	cndbicg_m.o ccfspt_m.o \
	cndbicg_mc.o ccfspt_mc.o 

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
