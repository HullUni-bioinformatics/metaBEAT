# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update
RUN apt-get update -y -qq --fix-missing
#RUN apt-get upgrade -y

RUN mkdir /home/metabarcoding
ENV PATH /home/metabarcoding:$PATH

RUN apt-get install ncbi-blast+
ADD external_software/pplacer /home/metabarcoding/
ADD external_software/guppy /home/metabarcoding/
ADD external_software/raxmlHPC-PTHREADS /home/metabarcoding/
ADD external_software/taxit /home/metabarcoding/
ADD external_software/hmmalign /home/metabarcoding/
ADD external_software/hmmbuild /home/metabarcoding/
