CXX = g++
RM = rm -f
CFLAGS = -pthread -std=c++0x -I./ -I./gzstream -m64 -Wall -g -Wno-sign-compare
DEBUGFLAGS   = -O0 -D _DEBUG
RELEASEFLAGS = -O2 -D NDEBUG 
TARGET  = optimap

SOURCES=main.cpp indexing.cpp gzstream.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(RELEASEFLAGS) -o $(TARGET) $(OBJECTS) -lz

main.o: main.cpp constants.h indexing.h types.h logger.h string_functions.h stats.h common.h
	$(CXX) $(CFLAGS) $(RELEASEFLAGS) -c main.cpp

indexing.o: indexing.cpp indexing.h
	$(CXX) $(CFLAGS) $(RELEASEFLAGS) -c indexing.cpp

gzstream.o: gzstream/gzstream.cpp gzstream/gzstream.h
	$(CXX) $(CFLAGS) $(RELEASEFLAGS) -c gzstream/gzstream.cpp

clean:
	$(RM) $(OBJECTS)

dist-clean:
	$(RM) $(OBJECTS) $(TARGET)
