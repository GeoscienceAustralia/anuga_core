###
#
# Config

CC = gcc
CFLAGS = 
OBJS = random.o
AR = ar rv
RANLIB = ranlib

## Merge the extra functions specified here into the metis library
merge: $(OBJS) ../$(METIS_DIR)/libmetis.a
	$(AR) ../$(METIS_DIR)/libmetis.a $(OBJS)
	$(RANLIB) ../$(METIS_DIR)/libmetis.a

../$(METIS_DIR)/libmetis.a:
	make -C ../$(METIS_DIR)/Lib CC=gcc

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJS)
	-rm -f *~
