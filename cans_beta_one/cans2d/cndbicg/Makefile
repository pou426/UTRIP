CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = gtnorm.o bicgstab1.o bicgstab2.o \
	iludcmp.o psolv.o cndprm.o \
	cndbicg.o ccfspt.o ccfunf.o \
	cndbicg_m.o ccfspt_m.o cndbicg_m3.o ccfspt_m3.o \
	cndbicg_c.o ccfspt_c.o ccfunf_c.o \
	cndbicg_mc.o ccfspt_mc.o cndbicg_m3c.o ccfspt_m3c.o \
	cndbicg_s.o ccfspt_s.o ccfunf_s.o \
	cndbicg_ms.o ccfspt_ms.o cndbicg_m3s.o ccfspt_m3s.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
