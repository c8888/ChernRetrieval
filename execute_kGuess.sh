#!/bin/sh

#PBS -N chern_kGuess

#PBS -q vlong

#PBS -l cput=240:00:00

#PBS -m abe

cd $HOME/ChernRetrieval/
/usr/local/math10/Executables/math -script kGuess.m
