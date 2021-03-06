ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
#INCDIR = -I./include -I$(shell echo ${G4INCLUDE})
vpath %.h ./include
vpath %.o obj
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) #$(INCDIR)

all: anaKK anaTrack AnaCuts AnaCutflow EffAndCross BackgroundAna anadimubck ShowBABAR\
     Ana_2D Ana_truthang Ana_OptEp EPoptimize2 MCshapeFit Ana_openAng ComparePolarAng2 Ana_simdata Ana_simdata_roAxis\
     CompareOpenAng_smear

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
