workspace "daisy"
configurations {"debug", "release"}

local mcu = "-mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard"

project "daisy_bitstream_autocorrelation"
kind "ConsoleApp"
toolset "gcc"
language "C++"
cppdialect "C++14"
targetdir "build/bin/%{cfg.buildcfg}"
objdir "build/obj/%{cfg.buildcfg}"
includedirs {"dependencies/libDaisy", "dependencies/libDaisy/src", "dependencies/libDaisy/src/sys",
             "dependencies/libDaisy/src/usbd", "dependencies/libDaisy/Drivers/CMSIS/Include",
             "dependencies/libDaisy/Drivers/CMSIS/DSP/Include",
             "dependencies/libDaisy/Drivers/CMSIS/Device/ST/STM32H7xx/Include",
             "dependencies/libDaisy/Drivers/STM32H7xx_HAL_Driver/Inc",
             "dependencies/libDaisy/Middlewares/ST/STM32_USB_Device_Library/Core/Inc",
             "dependencies/libDaisy/Middlewares/Third_Party/FatFs/src"}

gccprefix "arm-none-eabi-"
links {"daisy", "c", "m", "nosys"}
libdirs {"dependencies/libDaisy/build"}

files {"../src/**.cpp"}

buildoptions {mcu, "-Wall", "-fno-exceptions", "-fasm", "-fno-rtti", "-finline", "-finline-functions-called-once",
              "-fshort-enums", "-fno-move-loop-invariants", "-fno-unwind-tables", "-Wno-register"}

linkoptions {mcu, "--specs=nano.specs", "--specs=nosys.specs", "-Wl,--gc-sections", "-u _printf_float",
             "-Tdependencies/libdaisy/core/STM32H750IB_flash.lds"}

filter "configurations:debug"
defines {"DEBUG"}
symbols "On"

filter "configurations:release"
defines {"RELEASE"}
defines {"NDEBUG"}
optimize "On"
