CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = uwd_b.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
