# Makefile for bwt

# sources
OBJS				= divsufsort.o 
TARGET				= libdivsufsort.a
MAKEFILE			= Makefile

# options
CC					= gcc 
CFLAGS  	 		= -g -O3 -fomit-frame-pointer -fstrict-aliasing 
CPPFLAGS			= -Wall 
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=
AR       = /usr/bin/ar rc
RANLIB   = ranlib


# beautification magic stolen from git (added within hmmer source)
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_AR      = @echo '    ' AR $@;
endif


# targets
all: $(TARGET) 
$(TARGET): $(OBJS)

.c.o:
	${QUIET_CC}${CC} ${CFLAGS} ${CPPFLAGS} -o $@ -c $<


libdivsufsort.a: $(OBJS)
	${QUIET_AR}${AR} libdivsufsort.a $(OBJS)
	@${RANLIB} libdivsufsort.a
	@chmod 644 libdivsufsort.a


clean:
	$(RM) $(TARGET) $(OBJS)

# dependencies
$(OBJS): divsufsort.h $(MAKEFILE)
