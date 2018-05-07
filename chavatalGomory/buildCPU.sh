#! /bin/sh
#
# buildLAHCOpt.sh
# Copyright (C) 2016 haroldo <haroldo@soyuz>
#
# Distributed under terms of the MIT license.
#


CC=gcc
CPP=g++


CFLAGS="-Ofast -flto -fopenmp -std=c99 -I/opt/gurobi752/linux64/include/"
CPPFLAGS="-Ofast -flto -fopenmp -DGRB -I/opt/gurobi752/linux64/include/ "
LDFLAGS="-Ofast -flto -fopenmp -lpthread -lm -L/opt/gurobi752/linux64/lib/ -lgurobi75"
CSOURCES="cut_gpu.c Instance.c prepareCPU.c solutionGpu.c "
CPPSOURCES=" lp.cpp main.cpp"
BINDIR=./bin/opt/

rm $BINDIR/*


lnkFiles=""
echo building C sources ...
for cs in $CSOURCES;
do
    command="${CC} ${CFLAGS} -c $cs -o ${BINDIR}/`basename $cs .c`.o"
    printf "\t$command\n"
    lnkFiles="${lnkFiles}${BINDIR}/`basename $cs .c`.o "
    $command
done

echo building C++ sources ...
for cs in $CPPSOURCES;
do
    command="${CPP} ${CPPFLAGS} -c $cs -o ${BINDIR}/`basename $cs .cpp`.o"
    printf "\t$command\n"
    lnkFiles="${lnkFiles}${BINDIR}/`basename $cs .cpp`.o "
    $command
done

echo linking ...
    command="${CPP}  ${lnkFiles} ${LDFLAGS} -o ${BINDIR}/chavExecCPU"
    printf "\t$command\n"
    $command

cp bin/opt/chavExecCPU ./
