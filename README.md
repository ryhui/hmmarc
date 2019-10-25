# hmmarc
A hidden Markov model (HMM) to detect archaic segments in modern human genomes. The L-BFGS-B algorithm uses [Stephen Becker's implementation in C](https://github.com/stephenbeckr/L-BFGS-B-C).   
## Installation
```
git clone https://github.com/ryhui/hmmarc.git
cd hmmarc
make
```
## Usage
(Run `./hmmarc -h` for a summary of the options)
```
./hmmarc -m [model file] vit|post|est|bfgs [data file(s)]
```
**Model file**:

Includes parameters in the HMM (see example/model.hmm for an example).

N: number of hidden states;  
pi: starting state probabilities;  
T: time (in generations) since the admixture, used to calculate transition probabilities;  
n: number of observed states;  
emit: emission probabilities;  

In the decoding modes (vit/post), the model is used for inference; in the training modes (est/bfgs), the model specifies the initial parameter values.

**Data file**:

Records the observed allele-sharing states (see example/obs.HGDP00780.chr1.0.txt for an example).  
The file has three columns:

variant position;  
observed state;  
genetic distance (cM) from the previous variant.

**Mode**:

vit: output Viterbi sequence of the second state;  
est: output posterior probabilities of the second state at each input variant;  
est: Baum-Welch training;  
bfgs: numerical optimization using the L-BFGS-B algorithm.

## Examples
To print the Viterbi sequence of archaic segments:
```
./hmmarc -m example/model.hmm vit example/obs.HGDP00780.chr1.0.txt
```
To print the posterior probabilities of the archaic state at each input site:
```
./hmmarc -m example/model.hmm post example/obs.HGDP00780.chr1.0.txt
```
