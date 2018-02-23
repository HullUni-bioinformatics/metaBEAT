# building on Szitenberg/Reprophylo
FROM szitenberg/reprophylo

MAINTAINER Christoph Hahn <christoph.hahn@uni-graz.at>

# Make sure apt is up to date
RUN apt-get update -y --fix-missing && \
apt-get install -y openjdk-7-jre-headless wget ncbi-blast+ unzip
#RUN apt-get upgrade -y

RUN pip install biom-format==2.1.4 biopython==1.70

RUN mkdir /src
WORKDIR /src

#set up taxtastic v0.8.5
RUN git clone https://github.com/fhcrc/taxtastic.git && \
pip install -U pip && \
./taxtastic/dev/install_pysqlite.sh
RUN pip install taxtastic==0.8.5 && \
rm -rf pysqlite-2.8.3 pysqlite-2.8.3.tar.gz taxtastic

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

#Add binaries
ADD external_software/pplacer \
external_software/guppy \
external_software/rppr \
external_software/raxmlHPC-PTHREADS \
external_software/hmmalign \
external_software/hmmbuild \
external_software/flash \
external_software/vsearch \
external_software/trimmomatic-0.32.jar \
external_software/fastq_to_fasta \
external_software/fastx_clipper \
external_software/fastx_reverse_complement \
external_software/process_shortreads \
external_software/fastq_quality_trimmer \
external_software/fastq_quality_filter \
scripts/fetch_from_db.py \
scripts/metaBEAT_global.py \
scripts/jplace_to_biom.py /usr/bin/
#add the functions to a PYTHONPATH location
ADD scripts/metaBEAT_global_misc_functions.py /home/reprophylo/

#create taxonomy database and put in expected place
RUN taxit new_database ./taxonomy.db && \
rm ./taxdmp.zip && \
mv ./taxonomy.db /usr/bin/

RUN mkdir /home/working

#create mount point at /working
#RUN mkdir /home/working/IN
#RUN mkdir /home/working/OUT
VOLUME /home/working
#make /working the working directory in the image
WORKDIR /home/working
RUN chmod -R a+rw /home/working

