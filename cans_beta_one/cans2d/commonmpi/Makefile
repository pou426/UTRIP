CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = ix2igx.o p_bbtoaa.o p_bbtoaa_c.o p_bbtoaa_s.o p_bbtoaa_t.o \
       exc_9.o exc_8.o exc_7.o exc_6.o exc_4.o exc_3.o exc_1.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
