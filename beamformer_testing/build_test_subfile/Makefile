CC = gcc
CFLAGS = -O2 -Wall -Wno-maybe-uninitialized
LIBS = -lfftw3 -lm

TARGETS = make_delayed_tones_plus_noise make_test_subfile
SOURCE1 = make_delayed_tones_plus_noise.c
SOURCE2 = make_test_subfile.c

all: $(TARGETS)

make_delayed_tones_plus_noise: $(SOURCE1)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

make_test_subfile: $(SOURCE2)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -f $(TARGETS)
