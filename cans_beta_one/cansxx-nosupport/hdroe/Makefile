CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = roe_mh.o roereset.o roeflux2_m.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
