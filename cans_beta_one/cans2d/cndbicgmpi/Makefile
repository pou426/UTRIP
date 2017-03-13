CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = p_bicgstab2.o p_cndbicg.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
