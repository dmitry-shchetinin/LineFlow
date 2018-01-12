CC = gcc
CFLAGS = -I./include -fPIC -O2 -g -Wall -Wno-unused-variable

SOURCES_LIB = $(shell echo src/line_flow.c)
OBJECTS_LIB = $(SOURCES_LIB:.c=.o)
TARGET_LIB = $(PWD)/lib/Unix/libline_flow.so

SOURCES_TEST = $(shell echo tests/*.c)
OBJECTS_TEST = $(SOURCES_TEST:.c=.o)
TARGET_TEST = $(SOURCES_TEST:.c=.out)

ifeq ($(OS),Windows_NT)
  OS_DETECTED := Windows
else
  OS_DETECTED := $(shell uname -s)
endif
$(warning $(OS_DETECTED))

ifeq ($(OS_DETECTED),Darwin)
	LDFLAGS +=-dynamiclib -Wl,-install_name,$(TARGET_LIB)
else
	LDFLAGS += -shared
endif

.PHONY: all
all : lib

.PHONY: lib
lib : $(TARGET_LIB)

$(TARGET_LIB) : $(OBJECTS_LIB)
	mkdir -p ./lib
	$(CC) $(CFLAGS) -o $@ $(OBJECTS_LIB) $(LDFLAGS) -lm

.PHONY: test
test : $(TARGET_TEST)
tests/%.out: tests/%.c
	$(CC) $(CFLAGS) -I./include -L./lib/Unix -Wl,-rpath ./lib/Unix -o $@ $< -lline_flow -lm
	./tests/test_line_flow.out
	rm ./tests/test_line_flow.out

.PHONY: clean
clean :
	rm -f $(OBJECTS_LIB) $(TARGET_LIB)
	rm -f $(OBJECTS_TEST) $(TARGET_TEST)
