# Makefile for bwt

# sources
OBJS					= divsufsort.o 
TARGET				= libbwt.a
MAKEFILE			= Makefile

# options
CC						= gcc
#OUTPUT_OPTION	= -o $@
#CFLAGS				= -O3 -fomit-frame-pointer
CFLAGS				= -g
CPPFLAGS			= -Wall #-DNDEBUG
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=
AR       = /usr/bin/ar rcv
RANLIB   = ranlib


# targets
all: $(TARGET) 
$(TARGET): $(OBJS)


libbwt.a: $(OBJS)
	$(AR) libbwt.a $(OBJS)
	$(RANLIB) libbwt.a
	chmod 644 libbwt.a


distclean: clean
clean:
	$(RM) $(TARGET) $(OBJS)

# dependencies
$(OBJS): divsufsort.h $(MAKEFILE)
