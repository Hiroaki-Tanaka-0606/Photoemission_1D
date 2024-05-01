# cpp builder
FROM ubuntu:latest AS builder
RUN apt-get -y update
RUN apt-get -y install gfortran emacs build-essential wget
WORKDIR /source/lapack
RUN wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.0.tar.gz && tar xzf v3.12.0.tar.gz
WORKDIR /source/lapack/lapack-3.12.0/
RUN cp make.inc.example make.inc && make

WORKDIR /source/cpp
COPY ./cpp /source/cpp
RUN make

FROM busybox:latest
COPY --from=builder /source/cpp/*.o /work/
WORKDIR /work
CMD ["./photoemission_well.o"]