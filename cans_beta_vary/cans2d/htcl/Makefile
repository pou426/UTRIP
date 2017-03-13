CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = htclprm.o htcl.o cooldef.o htstdef.o htclad.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
