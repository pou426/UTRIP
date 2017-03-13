CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = \
cip_m_bg.o cip_mt_bg.o ciprdy_m_b.o \
cip_sh.o  cippost_sh.o  cippre_sh.o  ciprdy_sh.o  moc_sm.o \
cip_sm.o  cippost_sm.o  cippre_sm.o  ciprdy_sm.o  moclag_sm.o


'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
