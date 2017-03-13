FC = mpif77
CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = p_cndsor.o   p_cndsor_m.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
