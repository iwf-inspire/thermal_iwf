################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../benchmarks/ex1_ex2_reconstruction.c \
../benchmarks/ex3_transient_heat.c \
../benchmarks/ex4_steady_heat.c \
../benchmarks/ex5_eng_model.c 

OBJS += \
./benchmarks/ex1_ex2_reconstruction.o \
./benchmarks/ex3_transient_heat.o \
./benchmarks/ex4_steady_heat.o \
./benchmarks/ex5_eng_model.o 

C_DEPS += \
./benchmarks/ex1_ex2_reconstruction.d \
./benchmarks/ex3_transient_heat.d \
./benchmarks/ex4_steady_heat.d \
./benchmarks/ex5_eng_model.d 


# Each subdirectory must supply rules for building sources it contributes
benchmarks/%.o: ../benchmarks/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=c99 -std=gnu99 -DNDEBUG -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


