FROM ubuntu:22.04

# disable prompt during packages installation
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y make \
                   build-essential \
                   git \
                   cmake \
                   libboost-all-dev 
# default-jdk if a baseline in Java is needed
# vim gdb for debug

# remove package caches
RUN rm -rf /var/lib/apt/lists/*
RUN apt clean
