# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update
RUN apt-get update -y -qq --fix-missing
#RUN apt-get upgrade -y

RUN mkdir /home/metabarcoding
#RUN mkdir /home/bin
#ENV PATH /home/bin:$PATH

RUN apt-get install ncbi-blast+
ADD external_software/pplacer /usr/bin/
ADD external_software/guppy /usr/bin/
ADD external_software/raxmlHPC-PTHREADS /usr/bin/
ADD external_software/taxit /usr/bin/
ADD external_software/hmmalign /usr/bin/
ADD external_software/hmmbuild /usr/bin/
ADD external_software/flash /usr/bin/
ADD external_software/vsearch /usr/bin/
ADD external_software/trimmomatic-0.32.jar /usr/bin/

#create mount point at /home
VOLUME /home	
#make /home the working directory in the image
WORKDIR /home	
RUN chmod -R a+rw /home

