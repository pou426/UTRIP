CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = cip_m.o    cipdxsrc.o	ctranspty.o  moclagx.o	mocx.o \
    cip_m_g.o  ciprdy_m.o	ctransptz.o  moclagy.o	mocy.o \
    cipadv.o   ctransptx.o	intpvl.o     moclagz.o	mocz.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS) 
