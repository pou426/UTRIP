CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = tvdminmod.o \
     roeflux_a.o roe_a.o \
     roeflux_h.o roe_h.o roe_h_g.o roe_h_c.o roe_h_cg.o \
     roeflux_ht.o roe_ht.o roe_ht_g.o roe_ht_c.o roe_ht_cg.o \
     roeflux_m.o roe_m.o roe_m_g.o roe_m_bg.o \
     roeflux_m3.o roe_m3.o roe_m3_g.o \
     roeflux_mt.o roe_mt.o roe_mt_g.o roe_mt_bg.o \
     roeflux_m3t.o roe_m3t.o roe_m3t_g.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
