SRCS = $(shell find . -name "*.cpp") 

OBJS = $(addsuffix .o, $(basename $(SRCS)))

EXEC = painless

LIBS = -lglucose -L../glucose/parallel \
       -lminisat -L../minisat/build/release/lib/ \
       -llgl -L../lingeling/ \
       -lmapleCOMSPS -L../mapleCOMSPS/build/release/lib/ \
	   -lm4ri -L../mapleCOMSPS/m4ri-20140914/.libs \
       -lmaplesat -L../maplesat/build/release/lib/ \
       -lpthread -lz -lm #-static

CXXFLAGS = -I../glucose -I../minisat -I../lingeling -I../mapleCOMSPS -I../maplesat \
         -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -std=c++11 -O3

release: CXXFLAGS += -D NDEBUG
release: $(EXEC)

debug: CXXFLAGS += -DEBUG -g
debug: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC)
