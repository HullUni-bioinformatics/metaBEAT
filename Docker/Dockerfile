# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update -y --fix-missing
RUN apt-get install -y openjdk-7-jre-headless
#RUN apt-get upgrade -y

#RUN mkdir /home/metabarcoding
RUN pip install taxtastic==0.5.3
RUN pip install biom-format==2.1.4

RUN apt-get install ncbi-blast+
ADD external_software/pplacer /usr/bin/
ADD external_software/guppy /usr/bin/
ADD external_software/rppr /usr/bin/
ADD external_software/raxmlHPC-PTHREADS /usr/bin/
ADD external_software/hmmalign /usr/bin/
ADD external_software/hmmbuild /usr/bin/
ADD external_software/flash /usr/bin/
ADD external_software/vsearch /usr/bin/
ADD external_software/trimmomatic-0.32.jar /usr/bin/
ADD external_software/fastq_to_fasta /usr/bin/ 
ADD external_software/fastx_clipper /usr/bin/
ADD external_software/fastx_reverse_complement /usr/bin/
ADD external_software/process_shortreads /usr/bin/

RUN mkdir /home/working
RUN mkdir /src
ENV PATH /src:$PATH
RUN taxit new_database -d /src/taxonomy.db
RUN rm /src/taxdmp.zip
RUN mv /src/taxonomy.db /usr/bin/

ADD scripts/fetch_from_db.py /usr/bin/
ADD scripts/metaBEAT.py /usr/bin/
ADD scripts/jplace_to_biom.py /usr/bin/

ADD external_software/fastq_quality_trimmer /usr/bin/
ADD external_software/fastq_quality_filter /usr/bin/

#create mount point at /working
RUN mkdir /home/working/IN
RUN mkdir /home/working/OUT
VOLUME /home/working
#make /working the working directory in the image
WORKDIR /home/working	
RUN chmod -R a+rw /home/working

