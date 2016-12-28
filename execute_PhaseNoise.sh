#!/bin/sh

#PBS -N cn_PhNoise

#PBS -q mp16

#PBS -l cput=240:00:00

#PBS -m abe

cd $HOME/ChernRetrieval/
math -script cnAtPhaseNoise.m
