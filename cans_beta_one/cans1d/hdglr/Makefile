CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = monoto.o glr_h.o  glr_h_c.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
