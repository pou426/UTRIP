CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = mlwfull.o mlwhalf.o mlwartv.o mlwsrch.o mlwsrcf.o \
    mlw_a.o mlw_b.o \
    mlw_h.o mlw_h_g.o mlw_h_c.o mlw_h_cg.o \
    mlw_ht.o mlw_ht_g.o mlw_ht_c.o mlw_ht_cg.o \
    mlw_m.o mlw_m_g.o mlw_m_bg.o \
    mlw_m3.o mlw_m3_g.o \
    mlw_mt.o mlw_mt_g.o mlw_mt_bg.o \
    mlw_m3t.o mlw_m3t_g.o \
    mlw_em.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
