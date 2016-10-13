################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../intlib.f90 \
../landau-teller.f90 \
../landaumod.f90 

OBJS += \
./intlib.o \
./landau-teller.o \
./landaumod.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

intlib.o: ../intlib.f90

landau-teller.o: ../landau-teller.f90 intlib.o landaumod.o

landaumod.o: ../landaumod.f90


