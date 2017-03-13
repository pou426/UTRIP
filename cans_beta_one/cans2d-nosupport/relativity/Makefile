CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a

OBJS = rtnewt_sh.o rtnewt_sm.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
