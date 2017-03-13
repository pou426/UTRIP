CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = cipadv.o cipdxsrc.o intpvl.o moc.o moclag.o ctranspt.o \
	ciprdy_a.o cip_a.o \
	ciprdy_h.o cip_h.o cip_h_c.o cip_h_g.o cip_h_cg.o \
	ciprdy_ht.o cip_ht.o cip_ht_c.o cip_ht_g.o cip_ht_cg.o \
	ciprdy_m.o cip_m.o cip_m_g.o \
	ciprdy_m_b.o cip_m_bg.o \
	ciprdy_m3.o cip_m3.o cip_m3_g.o \
	ciprdy_mt.o cip_mt.o cip_mt_g.o \
	ciprdy_m3t.o cip_m3t.o cip_m3t_g.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
