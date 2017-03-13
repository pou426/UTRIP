CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a
OBJS = bdsmnx.o bdsmpx.o bdspnx.o bdsppx.o bdfrex.o bdcnsx.o bdperx.o \
	bmsppx.o bmsmpx.o bdfrdx.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
