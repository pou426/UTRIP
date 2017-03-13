CANS_LIB_DIR = ../..

TARGET = $(CANS_LIB_DIR)/libcans2d.a
OBJS = fill0.o melloc.o  relax.o   selfginit.o \
         addint.o  interp.o  mglin.o   resid.o   selfgper.o \
         copy.o    rstrct.o  slvsml.o ggtoggm.o

'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
