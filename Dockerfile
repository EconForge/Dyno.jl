FROM andrewosh/binder-base

# MAINTAINER Andrew Osheroff <andrewosh@gmail.com>

USER root

# Add Julia dependencies
RUN apt-get update && apt-get install -y wget && apt-get clean
RUN apt-get install -y curl libzmq3


USER main

RUN wget https://julialang.s3.amazonaws.com/bin/linux/x64/0.5/julia-0.5.0-linux-x86_64.tar.gz
RUN mkdir $HOME/julia
RUN tar xvf julia-0.5.0-linux-x86_64.tar.gz -C $HOME/julia --strip-components=1
ENV PATH $PATH:$HOME/julia/bin

# Install Julia kernel
RUN julia -e 'Pkg.add("IJulia")'
RUN julia -e 'Pkg.add("PyPlot")'

USER root
RUN touch /home/main/.curlrc
RUN chown main /home/main/.curlrc

USER main
RUN echo 'cacert=/etc/ssl/certs/ca-certificates.crt' > $HOME/.curlrc
RUN julia -e 'Pkg.add("Gadfly")'
