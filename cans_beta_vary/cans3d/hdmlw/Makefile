CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = mlwfull.o mlwhalf.o mlwartv.o mlwsrch.o mlwsrcf.o \
    mlw_h.o mlw_h_g.o mlw_h_c.o mlw_h_cg.o mlw_h_s.o mlw_h_sg.o \
    mlw_ht.o mlw_ht_g.o mlw_ht_c.o mlw_ht_cg.o mlw_ht_s.o mlw_ht_sg.o \
    mlw_m.o mlw_m_g.o mlw_m_e.o mlw_m_ge.o mlw_m_t.o \
    mlw_m_c.o mlw_m_cg.o mlw_m_ce.o mlw_m_cge.o \
    mlw_m_s.o mlw_m_sg.o mlw_m_se.o mlw_m_sge.o \
    mlw_mt.o mlw_mt_g.o mlw_mt_e.o mlw_mt_ge.o mlw_mt_t.o \
    mlw_mt_c.o mlw_mt_cg.o mlw_mt_ce.o mlw_mt_cge.o \
    mlw_mt_s.o mlw_mt_sg.o mlw_mt_se.o mlw_mt_sge.o \
    mlw_sh.o mlw_sm.o mlw_h_o.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
