ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
#INCDIR = -I./include -I$(shell echo ${G4INCLUDE})
vpath %.h ./include
vpath %.o obj
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) #$(INCDIR)

all: anaKK anaTrack AnaCuts AnaCutflow EffAndCross CrossWithTMinuit CrossWithTMinuit3 BackgroundAna anadimubck ShowBABAR\
     Ana_epadd Ana_epmin Ana_thetaadd Ana_thetamin Ana_tofadd Ana_tofmin Ana_2D Ana_truthang Ana_OptEp EPoptimize2

anaKK:Ana.C
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@
	@#echo "cleaning trash ..."
	@#-rm *.o


#$(OBJS): %.o: %.C
obj/%.o: %.C
	@echo "making object $@"
	@g++ $(ROOTINCLUDE) $(INCDIR) -c $< -o $@

%: %.C
	@echo "compiling $@"
	$(CC) `root-config --cflags --libs` -lRooFitCore -lRooFit -lMathMore -lMinuit $^ -o $@

.PHONY:clean
clean:
	-rm -f *.o src/*.o obj/*.o anaKK

cleanall:
	-rm -f anaKK *.o src/*.o *.eps *.pdf *.ps
