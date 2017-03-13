CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = cfl_a.o cfl_h.o cfl_ht.o cfl_m.o cfl_m_e.o cfl_m3.o cfl_m3_e.o \
       cfl_mt.o cfl_mt_e.o cfl_m3t.o cfl_m3t_e.o \
       cfl_h_s.o cfl_ht_s.o cfl_m_s.o cfl_m_se.o cfl_m3_s.o cfl_m3_se.o \
       cfl_mt_s.o cfl_mt_se.o cfl_m3t_s.o cfl_m3t_se.o \
       cfl_h_d.o \
	prtote.o tetopr.o \
	bbtocx.o bbtocy.o bbtocz.o \
	bbtocx_c.o bbtocy_c.o bbtocz_c.o \
	bbtocx_s.o bbtocy_s.o bbtocz_s.o \
	bbtoaa.o bbtoaa_c.o bbtoaa_s.o bbtoaa_t.o \
	grdrdy.o chkda.o chkdav.o dampx.o dampy.o saveini.o rangen.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
