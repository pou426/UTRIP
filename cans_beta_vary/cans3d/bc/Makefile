CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = bdsmnx.o bdsmpx.o bdspnx.o bdsppx.o bdfrex.o bdcnsx.o bdperx.o \
       bmsppx.o bmsmpx.o bdfrdx.o \
       bdsmny.o bdsmpy.o bdspny.o bdsppy.o bdfrey.o bdcnsy.o bdpery.o \
       bmsppy.o bmsmpy.o bdfrdy.o \
       bdsmnz.o bdsmpz.o bdspnz.o bdsppz.o bdfrez.o bdcnsz.o bdperz.o \
       bmsppz.o bmsmpz.o bdfrdz.o bdcndx.o bdcndy.o bdcndz.o bdiniz.o


'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
