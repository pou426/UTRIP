CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = lagpart.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
