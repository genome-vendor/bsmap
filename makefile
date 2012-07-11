CC=	g++

BIN = $(DESTDIR)/usr/bin
FLAGS= -DMAXHITS=1000 -DTHREAD -funroll-loops -Lsamtools -Isamtools -O3

OLIGOLEN= -DREAD_144
# options: -DREAD_48, -DREAD_80, -DREAD_144

LIBS=  
THREAD=	-lpthread

SOURCE = align dbseq main pairs param reads utilities
OBJS1= $(patsubst %,%.o,$(SOURCE))

all: bsmap
%.o:%.cpp
	$(CC) $(FLAGS) $(LIBS) $(REF_MODE) $(OLIGOLEN) -c $< -o $@
bsmap: $(OBJS1)
	(cd samtools; make)
	$(CC) $(FLAGS) $(LIBS) $(REF_MODE) $(OLIGOLEN) $^ -o $@ $(THREAD) -lbam -lz
	rm -f *.o

clean:
	rm -f *.o *~ bsmap
	(cd samtools; make clean)
install:
	install -d $(BIN)
	install ./bsmap $(BIN)
	install ./sam2bam.sh $(BIN)
	install ./methratio.py $(BIN)
