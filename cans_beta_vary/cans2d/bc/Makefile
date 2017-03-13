CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = bdsmnx.o bdsmpx.o bdspnx.o bdsppx.o bdfrex.o bdcnsx.o bdperx.o \
	bmsppx.o bmsmpx.o bdfrdx.o bdinix.o \
       bdsmny.o bdsmpy.o bdspny.o bdsppy.o bdfrey.o bdcnsy.o bdpery.o \
	bmsppy.o bmsmpy.o bdfrdy.o bdcndy.o bdiniy.o bdfdvx.o bdfdvy.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
