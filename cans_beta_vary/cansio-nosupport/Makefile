CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcansio.a
OBJS =  ncopn1d.o ncopn2d.o ncopn3d.o \
        ncgetattc.o ncputattc.o \
	ncdefss.o ncputss.o ncgetss.o \
	ncdefssi.o ncputssi.o ncgetssi.o \
	ncdefs1.o ncputs1.o ncgets1.o \
	ncdefs2.o ncputs2.o ncgets2.o \
	ncdefs3.o ncputs3.o ncgets3.o \
	ncputosc.o ncgetosc.o ncputosi.o ncgetosi.o ncputos.o ncgetos.o \
	ncputo1.o ncgeto1.o ncputo2.o ncgeto2.o ncputo3.o ncgeto3.o

'.update' : $(OBJS)
	ar -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)  ../libcansio.a
cleanall :
	make clean
