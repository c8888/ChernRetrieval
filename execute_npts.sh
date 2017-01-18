#!/bin/sh

#PBS -N cn_npts

#PBS -q vlong

#PBS -l cput=240:00:00

#PBS -m abe

cd $HOME/ChernRetrieval/
math -script cnAtnptsBZ.m
