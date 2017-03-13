CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = residue.o sorbr.o cndprm.o \
cndsor.o ccfspt.o ccfunf.o cndsor_c.o ccfspt_c.o ccfunf_c.o \
cndsor_m.o ccfspt_m.o cndsor_mc.o ccfspt_mc.o \
cndsor_m3.o ccfspt_m3.o cndsor_m3c.o ccfspt_m3c.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
