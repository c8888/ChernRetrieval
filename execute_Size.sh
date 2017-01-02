#!/bin/sh

#PBS -N cn_Size

#PBS -q mp8

#PBS -l cput=340:00:00

#PBS -m abe

cd $HOME/ChernRetrieval/
math -script cnAtSize.m
