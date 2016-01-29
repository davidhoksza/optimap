### Optimap - tool for optical mapping


#### Installation

### Requirements

- gcc 4.6 or newer
- [htslib](https://github.com/samtools/htslib): C-library for handling high-throughput sequencing data

### Download optimap

```
git clone https://github.com/davidhoksza/optimap.git
```


### Build

```
make
```

#### Usage

Command-line usage

Run:

```
./omap -r data/ref.map -o data/experiment.maps -k 10 -b 0 -e 3
```
  
to map first 4 (parameters *b* and *e* - optional) maps in *experiment.maps* to *ref.map*, obtain 10 highest scoring alignments and return them on the standard output. The input map files can be either plain text files or gzipped files.

Run

```
./omap -h
```

to get the description of all command-line parameters.
