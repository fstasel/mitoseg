FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y git ca-certificates build-essential cmake libopencv-dev libboost-dev libboost-program-options-dev libyaml-cpp-dev --no-install-recommends && \
    apt-get clean

WORKDIR /app
RUN git clone https://github.com/fstasel/mitoseg.git .
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc)

RUN chmod +x /app/docker/entrypoint.sh
ENTRYPOINT ["/app/docker/entrypoint.sh"]
