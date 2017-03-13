CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a

OBJS = roe_me_g.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
