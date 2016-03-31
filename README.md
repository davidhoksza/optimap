### Optimap - tool for optical mapping


#### Installation

### Requirements

- gcc 4.6 or newer

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
./optimap -r data/ref.map -o data/experiment.maps -k 10 -b 0 -e 3 --read-error-stddev 0.05 --cut-eff 0.8
```
  
to map first 4 (parameters *b* and *e* - optional) maps in *experiment.maps* to *ref.map*, obtain 10 highest scoring alignments and return them on the standard output. The input map files can be either plain text files or gzipped files (file name needs to end with .gz).

To algorithm seeks to identify alignment with biggest probability. If we had a reference optical map (R), got an optical map (O) from the same sequence and the machine which produces optical maps was flawless, we would expect to find exactly one-to-one mapping between R and O. However, this is not the case and the machine produces three types of error which are modelled by three types of input parameters:

* Sizing error, i.e. how inprecise the lengths of the fragments are. For a fragment of length R, this error is modelled as N(0, (*s*R)^2) . The s parameter of the distribution can be set using *read-error-stddev*
* Cut (digestion) efficiency, i.e. how often a restriction site is missed -> the probabily of missing N restriction sites is (1 - *cut-eff*)^N.
* False cuts rate,  i.e. how often a random cut appears. The probability of N cuts is modelled as Poisson distribution with mean = * *false-cut-p* x segment_length. This should be a very rare event and by default, *false-cut-p* is set to 0.00000001.
 
**Correct setting of the probability distributions parameters heavilly influences discovered alignments!**

Run

```
./omap -h
```

to get the description of all command-line parameters.
