# make file for Patterny

CC = gcc
CFLAGS = -O2 -fno-common -IfLPS2
LDFLAGS = -lm
SRCDIR = .
TARGETS = bandy pickbanding moduley blocky runny repeaty fLPS2/fLPS2 

BANDY_SOURCE = bandy.c
PICKBANDING_SOURCE = pickbanding.c
MODULEY_SOURCE = moduley.c
BLOCKY_SOURCE = blocky.c
RUNNY_SOURCE = runny.c
REPEATY_SOURCE  = repeaty.c
FLPS2_SOURCE = fLPS2/fLPS2.c

all: $(TARGETS) 

bandy: $(BANDY_SOURCE)
	$(CC) $(CFLAGS) $(BANDY_SOURCE) $(LDFLAGS) -o $@

pickbanding: $(PICKBANDING_SOURCE)
	$(CC) $(CFLAGS) $(PICKBANDING_SOURCE) $(LDFLAGS) -o $@

moduley: $(MODULEY_SOURCE)
	$(CC) $(CFLAGS) $(MODULEY_SOURCE) $(LDFLAGS) -o $@

blocky: $(BLOCKY_SOURCE)
	$(CC) $(CFLAGS) $(BLOCKY_SOURCE) $(LDFLAGS) -o $@

runny: $(RUNNY_SOURCE)
	$(CC) $(CFLAGS) $(RUNNY_SOURCE) $(LDFLAGS) -o $@

repeaty: $(REPEATY_SOURCE)
	$(CC) $(CFLAGS) $(REPEATY_SOURCE) $(LDFLAGS) -o $@

fLPS2/fLPS2: $(FLPS2_SOURCE)
	$(CC) $(CFLAGS) $(FLPS2_SOURCE) $(LDFLAGS) -o $@


# Clean rule
clean:
	rm -f $(TARGETS) 

.PHONY: all clean $(TARGETS)




