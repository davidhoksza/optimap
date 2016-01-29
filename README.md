### Optimap - tool for optical mapping


#### Installation

### Requirements

gcc 4.6 or newer

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
./omap -r data/ref.map -o data/experiment.maps -m mapping.out -k 10 -b 0 -e 3
```
  
to map first 4 (parameters *b* and *e* - optional) maps in *experiment.maps* to *ref.map*, obtain 10 highest scoring alignments and store the results in *mapping.out*.

Run

```
./omap -h
```

to get the description of all command-line parameters.
