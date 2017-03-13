CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = cfl_a.o cfl_b.o cfl_h.o cfl_ht.o cfl_m.o cfl_m3.o cfl_mt.o cfl_m3t.o \
        prtote.o tetopr.o \
	grdrdy.o scrdy.o damp.o chkdav.o chkda.o rangen.o \
	cfl_em.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
