CC     = nvcc 
CFLAGS = -O3 -lineinfo --ptxas-options=-v --maxrregcount=64 -arch=compute_60 -code=sm_61  -lineinfo --Wno-deprecated-declarations --compiler-options "-msse4.1" 
HOME = ../
all:   px_cuda_packed px_cuda px_cuda_texture px_cuda_fp

px:     px.c Makefile
	${CC} ${CFLAGS} -o px -L ${HOME}/src -L ${HOME}/lib/easel -I ${HOME}/src -I /usr/include/mpi -I ${HOME}/lib/easel px.c -lhmmer -leasel -lm -lmpi

px_serial:  px_serial.c Makefile
	${CC} ${CFLAGS} -o px_serial -L ${HOME}/src -L ${HOME}/lib/easel -I ${HOME}/src -I /usr/include/mpi -I ${HOME}/lib/easel px_serial.c -lhmmer -leasel -lm -lmpi


engine:  engine.c Makefile
	${CC} ${CFLAGS} -o engine -L ${HOME}/src -L ${HOME}/lib/easel -I ${HOME}/src -I /usr/include/mpi -I ${HOME}/lib/easel engine.c -lhmmer -leasel -lm -lmpi

px_cuda: px_cuda.cu Makefile
	${CC} ${CFLAGS} -o px_cuda  px_cuda.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi

px_cuda_texture: px_cuda_texture.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_texture  px_cuda_texture.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi

px_cuda_simple: px_cuda_simple.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_simple  px_cuda_simple.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi


px_cuda_int: px_cuda_int.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_int  px_cuda_int.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi

px_cuda_int2: px_cuda_int2.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_int2  px_cuda_int2.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi

px_cuda_fp: px_cuda_fp.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_fp  px_cuda_fp.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi

px_cuda_packed: px_cuda_packed.cu Makefile
	${CC} ${CFLAGS} -o px_cuda_packed  px_cuda_packed.cu -L ${HOME}/src -L /usr/lib/openmpi -L ${HOME}/lib/easel -I /usr/include/mpi -I ${HOME}/src -I ${HOME}/lib/easel -lhmmer -leasel -lm -lmpi
clean:
	-rm *.o *~
	-rm px px_serial engine px_cuda
