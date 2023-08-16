A pipeline for BCR repertoire for - Roche 454 BCR mRNA with Multiplexed Samples. that were produced in the same fashion as those in [giang et al. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3699344/)


Library preperation and sequencing method:

The sequences were amplified using specific primers of the C-region and non coding V region of the IGH chain.


Each sequence is marked with a sample barcode at the beginning. The subsequent sequence can be either in the forward orientation, moving 5' to 3' along the V(D)J reading frame, or in the reverse complement orientation, moving in the opposite direction.


Input files:

1. fastq file of single-read sequencing
2. primer files

To test the pipeline:

We recommend to test the pipeline using a small example from the original reads that can download using the fastq-dump command.


```bash
fastq-dump -X 50000 SRR765688
```

And upload directly to dolphinnext. 


Output files:

1. {sampleName}_collapse-unique.fastq
2. {sampleName}_atleast-2.fastq
3. log tab file for each steps
4. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

* Quality control

	1. FilterSeq lenght
	2. FilterSeq quality
	
* Sample barcode and primer identification

	3. MaskPrimer score
	4. MaskPrimer align
	5. MaskPrimer align
	
* Deduplication and filtering

	6. CollapseSeq	
	7. SplitSeq group
	

Primers used:

* [MIDs](https://bitbucket.org/kleinstein/presto/src/master/examples/Jiang2013/SRR765688_MIDs.fasta)



* [CPrimers](https://bitbucket.org/kleinstein/presto/src/master/examples/Jiang2013/SRX190717_CPrimers.fasta)

* [VPrimers](https://bitbucket.org/kleinstein/presto/src/master/examples/Jiang2013/SRX190717_VPrimers.fasta)

