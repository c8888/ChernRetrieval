#!/bin/sh

#PBS -N chern_1

#PBS -q long

#PBS -l cput=72:00:00

#PBS -m abe

cd $HOME/ChernRetrieval/
math -script run.m
