CXX = g++
RM = rm -f
CFLAGS = -pthread -std=c++11 -I./ -m64 -Wall -g -Wno-sign-compare
DEBUGFLAGS   = -O0 -D _DEBUG
RELEASEFLAGS = -O2 -D NDEBUG
TARGET  = omap

SOURCES=main.cpp indexing.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(DEBUGFLAGS) -o $(TARGET) $(OBJECTS)

main.o: main.cpp constants.h indexing.h types.h logger.h string_functions.h
	$(CXX) $(CFLAGS) $(DEBUGFLAGS) -c main.cpp

indexing.o: constants.h types.h
	$(CXX) $(CFLAGS) $(DEBUGFLAGS) -c indexing.cpp

clean:
	$(RM) $(OBJECTS)

dist-clean:
	$(RM) $(OBJECTS) $(TARGET)
