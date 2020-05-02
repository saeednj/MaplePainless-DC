FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install wget -y
RUN apt-get install unzip
RUN apt-get install build-essential -y
RUN apt-get install zlib1g-dev libm4ri-dev -y
RUN DEBIAN_FRONTEND=noninteractive apt install -y iproute2 cmake python python-pip build-essential gfortran wget curl
RUN pip install supervisor awscli
RUN apt-get install openmpi-bin openmpi-common libopenmpi-dev iputils-ping -y
ADD painless painless

RUN cd painless && make -j 4
ADD run.sh supervised-scripts/run.sh
#ADD make_combined_hostfile.py supervised-scripts/make_combined_hostfile.py
RUN chmod 755 supervised-scripts/run.sh
EXPOSE 22

CMD supervised-scripts/run.sh
