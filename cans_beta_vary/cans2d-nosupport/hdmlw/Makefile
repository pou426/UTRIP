CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = mlw_sh.o mlw_sm.o mlw_h_o.o mlw_m_o.o \
    mlw_h_d.o mlw_h_dg.o 

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
