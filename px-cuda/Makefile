CC     = nvcc
CFLAGS = -O3 -g -lineinfo --compiler-options "-msse4.1" -gencode=arch=compute_61,code=compute_61
HOME = ../../
all:   px px_serial engine px_cuda px_cuda_texture

px:     px.c Makefile
	${CC} ${CFLAGS} -o px -L ${HOME}hmmer/src -L ${HOME}easel -I ${HOME}hmmer/src -I /usr/include/mpi -I ${HOME}easel px.c -lhmmer -leasel -lm -lmpi

px_serial:  px_serial.c Makefile
	${CC} ${CFLAGS} -o px_serial -L ${HOME}hmmer/src -L ${HOME}easel -I ${HOME}hmmer/src -I /usr/include/mpi -I ${HOME}easel px_serial.c -lhmmer -leasel -lm -lmpi


engine:  engine.c Makefile
	${CC} ${CFLAGS} -o engine -L ${HOME}hmmer/src -L ${HOME}easel -I ${HOME}hmmer/src -I /usr/include/mpi -I ${HOME}easel engine.c -lhmmer -leasel -lm -lmpi

px_cuda: px_cuda.cu Makefile
	${CC} ${CFLAGS} -o px_cuda  px_cuda.cu -L ${HOME}hmmer/src -L /usr/lib/openmpi -L ${HOME}easel -I /usr/include/mpi -I ${HOME}hmmer/src -I ${HOME}easel -lhmmer -leasel -lm -lmpi

px_cuda_texture: px_cuda_texture.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_texture  px_cuda_texture.cu -L ${HOME}hmmer/src -L /usr/lib/openmpi -L ${HOME}easel -I /usr/include/mpi -I ${HOME}hmmer/src -I ${HOME}easel -lhmmer -leasel -lm -lmpi
clean:
	-rm *.o *~
	-rm px px_serial engine px_cuda
