CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = rtbis_sh.o rtnewt_sm.o rtnewt1_sm.o rtnewt_sh.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
