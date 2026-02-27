CC = g++
CXXFLAGS += -I../include
CXXFLAGS += -DCUP_BACKEND_QUASI_STATIC
CFLAGS = -Wall -I/usr/include/ -I/usr/include/eigen3 -L/usr/local/lib
CFLAGS += $(CXXFLAGS)
LIBS = -lgsl -lgslcblas -lm -larmadillo

# List of test programs
TARGETS = time_behavior threshold steady_state

# Define how to build each target
all: $(TARGETS)

time_behavior: time_behavior.cxx
	$(CC) $(CFLAGS) time_behavior.cxx -o time_behavior $(LIBS)

threshold: threshold.cxx
	$(CC) $(CFLAGS) threshold.cxx -o threshold $(LIBS)

steady_state: steady_state.cxx
	$(CC) $(CFLAGS) steady_state.cxx -o steady_state $(LIBS)

clean:
	rm -f $(TARGETS)
