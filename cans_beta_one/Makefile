all : 
	for dir in cansio cans1d cans2d cans3d ; do \
	   (cd $$dir; make); \
	done

allmpi : 
	for dir in cansio cans1d cans2d cans3d ; do \
	   (cd $$dir; make mpi); \
	done


clean :
	for dir in cansio cans1d cans2d cans3d ; do \
	   (cd $$dir; make clean); \
	   (cd $$dir; make cleanmpi); \
	done

cleanmpi :
	for dir in cansio cans1d cans2d cans3d ; do \
	   (cd $$dir; make cleanmpi); \
	done

plot : 
	for dir in cansio cans1d cans2d cans3d ; do \
	   (cd $$dir; make plot); \
	done
