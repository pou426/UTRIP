CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = cfl_a.o cfl_h.o cfl_h_c.o cfl_h_s.o \
	cfl_m.o cfl_m_e.o cfl_m_c.o cfl_m_ce.o cfl_m_s.o cfl_m_se.o \
	cfl_ht.o cfl_ht_c.o cfl_ht_s.o \
	cfl_mt.o cfl_mt_e.o cfl_mt_c.o cfl_mt_ce.o cfl_mt_s.o cfl_mt_se.o \
	prtote.o tetopr.o \
	bbtocx.o bbtocy.o bbtocz.o \
	bbtocx_c.o bbtocy_c.o bbtocz_c.o \
	bbtocx_s.o bbtocy_s.o bbtocz_s.o \
	grdrdy.o grdrdy2.o chkda.o chkdav.o rangen.o saveini.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
