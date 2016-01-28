### Optimap - tool for optical mapping


#### Installation

### Requirements

gcc 4.7 or newer

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
./omap -r data/ref.map -o data/experiment.maps -m mapping.out -k 10
```
  
to map maps in **experiment.maps** to *ref.map*, obtain 10 highest scoring alignments and store the results in *mapping.out*.

Run

```
./omap -h
```

to get the description of all command-line parameters.
