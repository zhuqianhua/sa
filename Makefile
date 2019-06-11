GXX = g++ 
SLD = 
DLD = -lpthread 
OBJ = saopt.o sasub.o sanw.o sa.o
VPATH = src

sa : ${OBJ}
	${GXX} -o sa ${OBJ} ${DLD} ${SLD}
.PHONY : clean
clean : 
	rm -f ${OBJ}

saopt.o:saopt.h
sasub.o:sasub.h
sanw.o:sanw.h
sa.o:sa.h

