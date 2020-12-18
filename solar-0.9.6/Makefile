CC=			g++
CFLAGS=		-g -O3 -Wall #-DLIH_DEBUG
OBJS= 		sBlast.o readFile.o sortBlast.o crossBlast.o flag_rep.o write_result.o cluster.o

.SUFFIXES:.c .cc .o

.c.o:sBlast.h crossBlast.h
		$(CC) -c $(CFLAGS) $< -o $@
.cc.o:sBlast.h crossBlast.h
		$(CC) -c $(CFLAGS) $< -o $@

all:solar
solar:$(OBJS)
		$(CC) $(CFLAGS) $(OBJS) -o $@
		
clean:
		rm -f *.o *.yy.c solar gmon.out *.aux *.ps *.dvi *.log *.ps.gz *.ps~ *.tar.gz
		rm -f *.cp *.fn *.ky solar.log *.pg *.toc *.tp *.vr *.cps *.pdf
