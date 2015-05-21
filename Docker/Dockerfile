# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update -y --fix-missing
RUN apt-get install -y openjdk-7-jre-headless
#RUN apt-get upgrade -y

#RUN mkdir /home/metabarcoding
RUN pip install taxtastic==0.5.3

RUN apt-get install ncbi-blast+
ADD external_software/pplacer /usr/bin/
ADD external_software/guppy /usr/bin/
ADD external_software/raxmlHPC-PTHREADS /usr/bin/
ADD external_software/hmmalign /usr/bin/
ADD external_software/hmmbuild /usr/bin/
ADD external_software/flash /usr/bin/
ADD external_software/vsearch /usr/bin/
ADD external_software/trimmomatic-0.32.jar /usr/bin/

RUN mkdir /src
ADD scripts/metaBEAT.py /src/
ENV PATH /src:$PATH

#create mount point at /working
RUN mkdir /home/working
VOLUME /home/working
#make /working the working directory in the image
WORKDIR /home/working	
RUN chmod -R a+rw /home/working

