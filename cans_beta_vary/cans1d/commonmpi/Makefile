CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a
OBJS = ix2igx.o exc_3.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
