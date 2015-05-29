__Docker context for metaBEAT - metaBarcoding and Environmental DNA analyses tool__

More details about the pipeline can be found on our Github [page](https://github.com/HullUni-bioinformatics/metaBEAT).

We build on the docker image for [Reprophylo](https://registry.hub.docker.com/u/szitenberg/reprophylo/) running Ubuntu 14.04.

CONTACT: <chrisi.hahni@gmail.com>

Start up the container by running, e.g. the following command:

```bash
sudo docker run -v ~/Dropbox/Github/Docker/test:/home/working -t -i chrishah/metabeat /bin/bash
```

or do the following to e.g. mount two directories. One contains the input data and will be mounted as write only. THe second directory is mounted as output directory. 

```bash
INDIR=/home/chrishah/Dropbox/Github/Docker/test/James/data/
OUTDIR=/home/chrishah/Dropbox/Github/Docker/test/James/output/
sudo docker run -v $INDIR:/home/working/IN/:ro -v $OUTDIR:/home/working/OUT/ -t -i chrishah/metabeat /bin/bash
```
