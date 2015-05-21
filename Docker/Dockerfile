# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update
RUN apt-get update -y -qq --fix-missing
#RUN apt-get upgrade -y

RUN mkdir /home/metabarcoding
RUN mkdir /home/bin
ENV PATH /home/bin:$PATH

RUN apt-get install ncbi-blast+
ADD external_software/pplacer /home/bin/
ADD external_software/guppy /home/bin/
ADD external_software/raxmlHPC-PTHREADS /home/bin/
ADD external_software/taxit /home/bin/
ADD external_software/hmmalign /home/bin/
ADD external_software/hmmbuild /home/bin/

#create mount point at /home
VOLUME /home	
#make /home the working directory in the image
WORKDIR /home	
RUN chmod -R a+rw /home

