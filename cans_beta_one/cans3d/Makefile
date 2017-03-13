CANS_LIB_DIR = ..
TARGET = $(CANS_LIB_DIR)/libcans3d.a
all :  $(TARGET)

MODULE = bc cndbicg cndsor common hdmlw hdroe htcl
MODULEMPI = cndbicgmpi cndsormpi commonmpi


$(TARGET) : 
	for dir in $(MODULE) ; do \
	   (cd $$dir; make); \
	done

cleanall : clean cleanmd cleanmpi

clean :
	for dir in $(MODULE) ; do \
	   (cd $$dir; make clean); \
	done
	rm -f $(TARGET)

cleanmd :
	for dir in md_* ; do \
	   (cd $$dir; \
	       echo $$dir; \
	    for sdir in *; do \
	       if [ -d $$sdir ]; then \
		echo $$sdir; \
	         (cd $$sdir; make clean); \
	       fi \
	    done \
	    ); \
	done


runmd :
	for dir in md_* ; do \
	   (cd $$dir; \
	       echo $$dir; \
	    for sdir in *; do \
	       if [ -d $$sdir ]; then \
		echo $$sdir; \
	         (cd $$sdir; make run); \
	       fi \
	    done \
	    ); \
	done

plot :
	for dir in md_* ; do \
	   (cd $$dir; \
	       echo $$dir; \
	    for sdir in *; do \
	       if [ -d $$sdir ]; then \
		echo $$sdir; \
	         (cd $$sdir; make plot); \
	       fi \
	    done \
	    ); \
	done
mpi :
	for dir in $(MODULEMPI) ; do \
	   (cd $$dir; make); \
	done

cleanmpi :
	for dir in $(MODULEMPI) ; do \
	   (cd $$dir; make clean); \
	done

