CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = tvdminmod.o tvdminmod1.o \
     roeflux_a.o roe_a.o \
     roeflux_h.o roe_h.o roe_h_g.o \
     roeflux_m.o roe_m.o roe_m_g.o roe_m_e.o \
     roeflux_m3.o roe_m3.o \
     roeflux_ht.o roe_ht.o roe_ht_g.o \
     roeflux_mt.o roe_mt.o \
     roeflux_m3t.o roe_m3t.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
