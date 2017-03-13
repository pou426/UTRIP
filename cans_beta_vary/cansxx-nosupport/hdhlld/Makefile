.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $<
CANS_LIB_DIR = ../..
TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS = hlld_flux.o hlld_mh.o


'.update' : $(OBJS)
	$(AR) -r $(TARGET) $(OBJS)
	touch .update
clean :
	rm -f $(OBJS)
