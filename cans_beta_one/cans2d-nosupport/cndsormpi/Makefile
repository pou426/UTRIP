CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = p_cndsor_m.o  p_cndsor_m3.o \
       p_cndsor_mc.o p_cndsor_m3c.o \
       p_cndsor_ms.o p_cndsor_m3s.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
