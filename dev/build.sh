#!/bin/bash
OPTION=$1

premake5 gmake2 || exit 1

if [ "$OPTION" = "clean" ]
then
    echo Cleaning project ...
    rm -fR build && rm Makefile && rm *.make || exit 1
    # make again after cleaning.
    premake5 gmake2 || exit 1
    shift
fi

#make verbose=1 -j7 || exit 1
make -j7 || exit 1

openocd -s /usr/local/share/openocd/scripts -f interface/stlink.cfg -f target/stm32h7x.cfg \
	-c "program ./build/bin/debug/daisy_bitstream_autocorrelation verify reset exit"