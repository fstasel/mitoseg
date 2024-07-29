FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y git ca-certificates build-essential cmake libopencv-dev libboost-dev libyaml-cpp-dev --no-install-recommends && \
    apt-get clean

WORKDIR /app
RUN git clone https://github.com/fstasel/mitoseg.git .
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j $(nproc)

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["entrypoint.sh"]
