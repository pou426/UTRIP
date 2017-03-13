CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans3d.a

OBJS = rtnewt_sm.o rtnewt_sh.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
