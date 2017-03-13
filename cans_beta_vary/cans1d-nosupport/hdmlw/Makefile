CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = mlw_m_o.o   mlw_sh.o    mlw_sh_o.o mlw_h_o.o    mlw_sm.o


'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
