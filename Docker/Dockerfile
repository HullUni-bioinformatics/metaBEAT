# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <chrisi.hahni@gmail.com>

# Make sure apt is up to date
RUN apt-get update -y --fix-missing
RUN apt-get install -y openjdk-7-jre-headless
RUN apt-get install -y wget
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
ADD external_software/fastq_quality_trimmer /usr/bin/
ADD external_software/fastq_quality_filter /usr/bin/

RUN mkdir /src
WORKDIR /src

#kraken
RUN wget -qO- wget https://github.com/DerrickWood/kraken/archive/v0.10.5-beta.tar.gz | tar xvz -C . && \
cd kraken-0.10.5-beta/ && \
./install_kraken.sh . && \
cd ..
ENV PATH /src/kraken-0.10.5-beta/:$PATH

#jellyfish
RUN wget -qO- !wget http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz | tar xvz -C . && \
cd jellyfish-1.1.11/ && \
./configure --prefix=$(pwd) && make && make install 
ENV PATH /src/jellyfish-1.1.11/bin:$PATH

#create taxonomy database
RUN taxit new_database -d /src/taxonomy.db && \
rm /src/taxdmp.zip && \
mv /src/taxonomy.db /usr/bin/

#add scripts
ADD scripts/fetch_from_db.py /usr/bin/
#ADD scripts/metaBEAT.py /usr/bin/
ADD scripts/metaBEAT_global.py /usr/bin/
ADD scripts/jplace_to_biom.py /usr/bin/
#add the functions to a PYTHONPATH location
ADD scripts/metaBEAT_global_misc_functions.py /home/reprophylo/

RUN mkdir /home/working

#create mount point at /working
#RUN mkdir /home/working/IN
#RUN mkdir /home/working/OUT
VOLUME /home/working
#make /working the working directory in the image
WORKDIR /home/working	
RUN chmod -R a+rw /home/working

