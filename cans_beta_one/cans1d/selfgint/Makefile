CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans1d.a

OBJS = getmc.o selfgint.o getmc_c.o selfgint_c.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
