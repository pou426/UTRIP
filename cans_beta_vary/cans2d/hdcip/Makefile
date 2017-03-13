CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a

OBJS = cipadv.o cipdxsrc.o intpvl.o moc.o moclag.o ctranspt.o \
        ctranspt_e.o cipb2cz.o \
	mocz.o moclagz.o ctransptz.o \
	cipbdvsym.o cipbdbsym.o cipbdv3sym.o cipbdb3sym.o \
	ciprdy_a.o cip_a.o \
	ciprdy_h.o cip_h.o cip_h_g.o \
	ciprdy_ht.o cip_ht.o cip_ht_g.o \
	ciprdy_m.o cip_m.o cip_m_g.o cip_m_e.o \
	ciprdy_mt.o cip_mt.o cip_mt_g.o \
	ciprdy_m3.o cip_m3.o cip_m3_g.o \
	ciprdy_m3t.o cip_m3t.o cip_m3t_g.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
