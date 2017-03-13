CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a

OBJS = monoto.o \
     glr_h_g.o \
     glr_h.o 

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
