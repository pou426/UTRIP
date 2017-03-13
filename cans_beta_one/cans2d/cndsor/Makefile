CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = gtnorm.o sorbr.o cndprm.o \
	cndsor.o ccfspt.o ccfunf.o \
	cndsor_m.o ccfspt_m.o cndsor_m3.o ccfspt_m3.o \
	cndsor_c.o ccfspt_c.o ccfunf_c.o \
	cndsor_mc.o ccfspt_mc.o cndsor_m3c.o ccfspt_m3c.o \
	cndsor_s.o ccfspt_s.o ccfunf_s.o \
	cndsor_ms.o ccfspt_ms.o cndsor_m3s.o ccfspt_m3s.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
