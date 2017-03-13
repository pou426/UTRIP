CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcansio.a
OBJS =  \
	dacdefparam.o dacdef0s.o dacdef1d.o dacdef1s.o \
	dacputparami.o dacputparamc.o dacputparamr.o dacputparamd.o \
	dacopnr1d.o dacopnr0s.o dacopnr1s.o dacopnr2s.o dacopnr3s.o \
	dacdef2d.o dacdef2s.o dacdef3d.o dacdef3s.o dacdef4s.o dacopnr4s.o \
	dacgetparamd.o dacgetparamr.o dacgetparami.o dacgetparamc.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)  ../libcansio.a
cleanall :
	make clean
