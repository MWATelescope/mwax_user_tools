# Handy recipies for DSPSR

## Installing DSPSR

Don't bother, it's too hard. Instead use my handy Docker file!

```shell
cd dspsr
docker build . -t dspsr
```

This will provide a ready to use dspsr suite which includes VDIF, Sigproc (filterbank), FITS and PSRDADA support.

## Handy quick start guide for DSPSR

Run the below commands inside the docker container using this:
```shell
$ docker run -it -v /data:/data dspsr:latest /bin/bash
psr@f0ce5259b795:~$ 
```
(The -v is to map a local volume to a directory in the container)

and the cd to where your data is located:
```shell
psr@f0ce5259b795:~$ cd /data/1456155320
psr@f0ce5259b795:/data/1456155320$ 
```

### Fold and coherently de-disperse

```shell
psr@f0ce5259b795:/data/1456155320$ dspsr -S 4 -D 40.938 -c 1.2381295438682 -b 128 -O 1456155320_ch109_beam01 1456155320_ch109_beam01.hdr
```
Where:
* -S -> skip the first N seconds (4 is used here due to quacktime)
* -D -> the known DM of the pulsar
* -c -> the known pulse period of the pulsar (seconds)
* -b -> number of phase bins to use
* -O -> output filename (.ar extension gets added by dspsr)
* FILENAME -> is the full path to the HDR file

### Stats

Execute the below to see lots of stats about your folded pulsar data:

```shell
psr@f0ce5259b795:/data/1456155320$ psrstat 1456155320_ch109_beam01.ar
```

### Plots

#### Pulse profile

```shell
psr@f0ce5259b795:/data/1456155320$ pav -DFTp 1456155320_ch109_beam01.ar -g 1456155320_ch109_beam01_profile.png/png
```

or

```shell
psr@f0ce5259b795:/data/1456155320$ psrplot -p flux -jFDp 1456155320_ch109_beam01.ar -D 456155320_ch109_beam01_profile2.png/png
```

#### Phase vs frequency (2D plot)
```shell
psr@f0ce5259b795:/data/1456155320$ psrplot -p freq -jTDp 1456155320_ch109_beam01.ar -D 456155320_ch109_beam01_phase_v_freq.png/png
```
NOTE: I could not get this to work- always get error 

#### Phase vs time
```shell
psr@f0ce5259b795:/data/1456155320$ psrplot -p time -jFDp 1456155320_ch109_beam01.ar -D 456155320_ch109_beam01_phase_v_time.png/png
```