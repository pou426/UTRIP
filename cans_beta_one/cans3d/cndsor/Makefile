CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = gtnorm.o sorbr.o cndprm.o \
	cndsor.o ccfspt.o cndsor_m.o ccfspt_m.o \
	cndsor_c.o ccfspt_c.o cndsor_mc.o ccfspt_mc.o 

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
