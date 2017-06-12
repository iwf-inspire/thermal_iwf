################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../methods/all_methods.c \
../methods/cspm.c \
../methods/nmfs.c \
../methods/pse.c \
../methods/rkpm.c \
../methods/sph.c 

OBJS += \
./methods/all_methods.o \
./methods/cspm.o \
./methods/nmfs.o \
./methods/pse.o \
./methods/rkpm.o \
./methods/sph.o 

C_DEPS += \
./methods/all_methods.d \
./methods/cspm.d \
./methods/nmfs.d \
./methods/pse.d \
./methods/rkpm.d \
./methods/sph.d 


# Each subdirectory must supply rules for building sources it contributes
methods/%.o: ../methods/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=c99 -std=gnu99 -DNDEBUG -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


