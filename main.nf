$HOSTNAME = ""
params.outdir = 'results'  


//* params.nproc =  1  //* @input @description:"How many processes to use for each step. Default 1"
//* params.edit_Filter_Sequence_Length_filter_seq =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"filter_seq"
//* params.edit_Filter_Sequence_Quality_filter_seq =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"filter_seq"
//* params.edit_MaskPrimer_1_maskprimer =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"MaskPrimers"
//* params.edit_MaskPrimer_2_maskprimer =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"MaskPrimers"
//* params.edit_MaskPrimer_3_maskprimer =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"MaskPrimers"
//* params.edit_collapse_sequences_collapse_seq =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"collapse_seq"
//* params.edit_split_seq_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"split_seq"
//* params.edit_Parse_header_table_parse_headers =  "no"  //* @dropdown @options:"yes","no"  @show_settings:"parse_headers"

//* autofill
if ($HOSTNAME == "default"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"

}

//* platform
if ($HOSTNAME == "ig03.lnx.biu.ac.il"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"
	$CPU  = 48
    $MEMORY = 300 
}
//* platform


//* autofill



if(params.edit_Filter_Sequence_Length_filter_seq == "no"){
    // Process Parameters for params.Filter_Sequence_Length_filter_seq_quality:
    params.Filter_Sequence_Length_filter_seq_quality.method = "length"
    params.Filter_Sequence_Length_filter_seq_quality.nproc = params.nproc
    params.Filter_Sequence_Length_filter_seq_quality.n_length = "300"
}

if(params.edit_Filter_Sequence_Quality_filter_seq == "no"){

    // Process Parameters for params.Filter_Sequence_Quality_filter_seq_quality:
    params.Filter_Sequence_Quality_filter_seq_quality.method = "quality"
    params.Filter_Sequence_Quality_filter_seq_quality.nproc = params.nproc
    params.Filter_Sequence_Quality_filter_seq_quality.q = "20"
}
if(params.edit_MaskPrimer_1_maskprimer == "no"){
	// Process Parameters for MaskPrimer_1_MaskPrimers:
	params.Mask_Primer_1_MaskPrimers.method = ["score"]
	params.Mask_Primer_1_MaskPrimers.start = [0]
	params.Mask_Primer_1_MaskPrimers.maxerror  = [0.1]
	params.Mask_Primer_1_MaskPrimers.mode  = ["cut"]
	params.Mask_Primer_1_MaskPrimers.primer_field = ["MID"] 
	params.Mask_Primer_1_MaskPrimers.nproc = params.nproc
	params.Mask_Primer_1_MaskPrimers.R1_primers = "${projectDir}/primers/SRR765688_MIDs.fasta"
}
if(params.edit_MaskPrimer_2_maskprimer == "no"){
	// Process Parameters for MaskPrimer_2_maskprimer:
	params.Mask_Primer_2_MaskPrimers.method = ["align"]
	params.Mask_Primer_2_MaskPrimers.maxlen = [50]
	params.Mask_Primer_2_MaskPrimers.maxerror  = [0.03]
	params.Mask_Primer_2_MaskPrimers.mode  = ["mask"]
	params.Mask_Primer_2_MaskPrimers.primer_field = ["VPRIMER"]
	params.Mask_Primer_2_MaskPrimers.nproc = params.nproc
	params.Mask_Primer_2_MaskPrimers.R1_primers = "${projectDir}/primers/SRX190717_VPrimers.fasta"
}
if(params.edit_MaskPrimer_3_maskprimer == "no"){
	// Process Parameters for MaskPrimer_3_maskprimer:
	params.Mask_Primer_3_MaskPrimers.method = "align"
	params.Mask_Primer_3_MaskPrimers.maxlen = 50
	params.Mask_Primer_3_MaskPrimers.maxerror  = 0.3
	params.Mask_Primer_3_MaskPrimers.mode  = "cut"
	params.Mask_Primer_3_MaskPrimers.primer_field = "CPRIMER"
	params.Mask_Primer_3_MaskPrimers.skiprc = "true"
	params.Mask_Primer_3_MaskPrimers.nproc = params.nproc
	params.Mask_Primer_3_MaskPrimers.R1_primers = "${projectDir}/primers/SRX190717_CPrimers.fasta"
}
if(params.edit_collapse_sequences_collapse_seq == "no"){
	// Process Parameters for collapse_sequences_collapse_seq:
	params.collapse_sequences_collapse_seq.max_missing = 20
	params.collapse_sequences_collapse_seq.inner = "true"
	params.collapse_sequences_collapse_seq.uf = "MID CPRIMER"
	params.collapse_sequences_collapse_seq.cf = "VPRIMER"
	params.collapse_sequences_collapse_seq.act = "set"
	params.collapse_sequences_collapse_seq.nproc = params.nproc
}

if(params.edit_split_seq_params == "no"){
    // Process Parameters for params.split_sequences_split_seq:
    params.split_sequences_split_seq.field = "DUPCOUNT"
    params.split_sequences_split_seq.num = 2
    params.split_sequences_split_seq.nproc = params.nproc
}
if(params.edit_Parse_header_table_parse_headers == "no"){
	// Process Parameters for Parse_header_table_parse_headers:
	params.Parse_header_table_parse_headers.method = "table"
	params.Parse_header_table_parse_headers.args = "ID DUPCOUNT MID CPRIMER VPRIMER"
	params.Parse_header_table_parse_headers.nproc = params.nproc
}


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_11_reads_g0_0}
 } else {  
	g_11_reads_g0_0 = Channel.empty()
 }

Channel.value(params.mate).into{g_12_mate_g0_0;g_12_mate_g0_11;g_12_mate_g0_7;g_12_mate_g2_9;g_12_mate_g2_12;g_12_mate_g2_11;g_12_mate_g3_9;g_12_mate_g3_12;g_12_mate_g3_11;g_12_mate_g4_9;g_12_mate_g4_12;g_12_mate_g4_11;g_12_mate_g1_7;g_12_mate_g1_5;g_12_mate_g1_0}


process Filter_Sequence_Length_filter_seq_quality {

input:
 set val(name),file(reads) from g_11_reads_g0_0
 val mate from g_12_mate_g0_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g0_0_reads0_g1_0
 set val(name), file("FS_*")  into g0_0_logFile1_g0_11
 set val(name), file("*_${method}-fail.fastq") optional true  into g0_0_reads22
 set val(name),file("out*") optional true  into g0_0_logFile3_g25_0

script:
method = params.Filter_Sequence_Length_filter_seq_quality.method
nproc = params.Filter_Sequence_Length_filter_seq_quality.nproc
q = params.Filter_Sequence_Length_filter_seq_quality.q
n_length = params.Filter_Sequence_Length_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Length_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed >> out_${R1}_FS.log
	"""
}


}


process Filter_Sequence_Length_parse_log_FL {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "FSL_log_tab/$filename"}
input:
 set val(name), file(log_file) from g0_0_logFile1_g0_11
 val mate from g_12_mate_g0_11

output:
 set val(name), file("*.tab")  into g0_11_logFile0_g0_7

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID LENGTH
"""
}


process Filter_Sequence_Length_report_filter_Seq_Lenght {

input:
 set val(name), file(log_files) from g0_11_logFile0_g0_7
 val matee from g_12_mate_g0_7

output:
 file "*.rmd"  into g0_7_rMarkdown0_g0_9


shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("lenght", 
	        paste("Minimun lenght for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the min lenght under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	lenght_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	lenght_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Lenght Scores
	
	Lenght filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq Lenght tool remove short reads with minimun lenght. 
	The most commonly used approach is to remove read with lenght below 35.
	
	```{r, echo=FALSE}
	plotFilterSeq(lenght_log_1, lenght_log_2, titles=plot_titles, cutoff=35 , sizing="figure")
	```
	
	`r figures("lenght")`
		
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("lenght", 
	        paste("Minimun lenght for",  plot_titles[1],
	              "The dotted line indicates the min lenght under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	lenght_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Lenght Scores
	
	Lenght filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq Lenght tool remove short reads with minimun lenght. 
	The most commonly used approach is to remove read with lenght below 35.
	
	```{r, echo=FALSE}
	plotFilterSeq(lenght_log_1, titles=plot_titles, cutoff=35 , sizing="figure")
	```
	
	`r figures("lenght")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Length_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "FSL_report/$filename"}
input:
 file rmk from g0_7_rMarkdown0_g0_9

output:
 file "*.html"  into g0_9_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g0_0_reads0_g1_0
 val mate from g_12_mate_g1_0

output:
 set val(name), file("*_${method}-pass.fastq")  into g1_0_reads00
 set val(name), file("FS_*")  into g1_0_logFile1_g1_5
 set val(name), file("*_${method}-fail.fastq") optional true  into g1_0_reads22
 set val(name),file("out*") optional true  into g1_0_logFile3_g25_0

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="quality"){
	q = "-q ${q}"
	n_length = ""
	n_missing = ""
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = ""
		n_length = ""
		n_missing = "-n ${n_missing}"
	}
}

readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed >> out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed >> out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed >> out_${R1}_FS.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "FSL_log_tab/$filename"}
input:
 set val(name), file(log_file) from g1_0_logFile1_g1_5
 val mate from g_12_mate_g1_5

output:
 set val(name), file("*.tab")  into g1_5_logFile0_g1_7

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

input:
 val matee from g_12_mate_g1_7
 set val(name), file(log_files) from g1_5_logFile0_g1_7

output:
 file "*.rmd"  into g1_7_rMarkdown0_g1_9


shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "FSQ_report/$filename"}
input:
 file rmk from g1_7_rMarkdown0_g1_9

output:
 file "*.html"  into g1_9_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_1_MaskPrimers {

input:
 val mate from g_12_mate_g2_11

output:
 set val(name), file("*_primers-pass.fastq")  into g2_11_reads0_g3_11
 set val(name), file("*_primers-fail.fastq") optional true  into g2_11_reads_failed11
 set val(name), file("MP_*")  into g2_11_logFile2_g2_9
 set val(name),file("out*")  into g2_11_logFile3_g25_0

script:
method = params.Mask_Primer_1_MaskPrimers.method
barcode_field = params.Mask_Primer_1_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_1_MaskPrimers.primer_field
barcode = params.Mask_Primer_1_MaskPrimers.barcode
revpr = params.Mask_Primer_1_MaskPrimers.revpr
mode = params.Mask_Primer_1_MaskPrimers.mode
failed = params.Mask_Primer_1_MaskPrimers.failed
nproc = params.Mask_Primer_1_MaskPrimers.nproc
maxerror = params.Mask_Primer_1_MaskPrimers.maxerror
umi_length = params.Mask_Primer_1_MaskPrimers.umi_length
start = params.Mask_Primer_1_MaskPrimers.start
extract_length = params.Mask_Primer_1_MaskPrimers.extract_length
maxlen = params.Mask_Primer_1_MaskPrimers.maxlen
skiprc = params.Mask_Primer_1_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_1_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_1_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process Mask_Primer_1_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP1_log_tab/$filename"}
input:
 val mate from g_12_mate_g2_9
 set val(name), file(log_file) from g2_11_logFile2_g2_9

output:
 set val(name), file("*.tab")  into g2_9_logFile0_g2_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_1_try_report_maskprimer {

input:
 set val(name), file(primers) from g2_9_logFile0_g2_12
 val matee from g_12_mate_g2_12

output:
 file "*.rmd"  into g2_12_rMarkdown0_g2_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_1_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "MP1_report/$filename"}
input:
 file rmk from g2_12_rMarkdown0_g2_16

output:
 file "*.html"  into g2_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_2_MaskPrimers {

input:
 val mate from g_12_mate_g3_11
 set val(name),file(reads) from g2_11_reads0_g3_11

output:
 set val(name), file("*_primers-pass.fastq")  into g3_11_reads0_g4_11
 set val(name), file("*_primers-fail.fastq") optional true  into g3_11_reads_failed11
 set val(name), file("MP_*")  into g3_11_logFile2_g3_9
 set val(name),file("out*")  into g3_11_logFile3_g25_0

script:
method = params.Mask_Primer_2_MaskPrimers.method
barcode_field = params.Mask_Primer_2_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_2_MaskPrimers.primer_field
barcode = params.Mask_Primer_2_MaskPrimers.barcode
revpr = params.Mask_Primer_2_MaskPrimers.revpr
mode = params.Mask_Primer_2_MaskPrimers.mode
failed = params.Mask_Primer_2_MaskPrimers.failed
nproc = params.Mask_Primer_2_MaskPrimers.nproc
maxerror = params.Mask_Primer_2_MaskPrimers.maxerror
umi_length = params.Mask_Primer_2_MaskPrimers.umi_length
start = params.Mask_Primer_2_MaskPrimers.start
extract_length = params.Mask_Primer_2_MaskPrimers.extract_length
maxlen = params.Mask_Primer_2_MaskPrimers.maxlen
skiprc = params.Mask_Primer_2_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_2_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_2_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process Mask_Primer_2_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP2_log_tab/$filename"}
input:
 val mate from g_12_mate_g3_9
 set val(name), file(log_file) from g3_11_logFile2_g3_9

output:
 set val(name), file("*.tab")  into g3_9_logFile0_g3_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_2_try_report_maskprimer {

input:
 set val(name), file(primers) from g3_9_logFile0_g3_12
 val matee from g_12_mate_g3_12

output:
 file "*.rmd"  into g3_12_rMarkdown0_g3_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_2_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "MP2_report/$filename"}
input:
 file rmk from g3_12_rMarkdown0_g3_16

output:
 file "*.html"  into g3_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Mask_Primer_3_MaskPrimers {

input:
 val mate from g_12_mate_g4_11
 set val(name),file(reads) from g3_11_reads0_g4_11

output:
 set val(name), file("*_primers-pass.fastq")  into g4_11_reads0_g5_16
 set val(name), file("*_primers-fail.fastq") optional true  into g4_11_reads_failed11
 set val(name), file("MP_*")  into g4_11_logFile2_g4_9
 set val(name),file("out*")  into g4_11_logFile3_g25_0

script:
method = params.Mask_Primer_3_MaskPrimers.method
barcode_field = params.Mask_Primer_3_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_3_MaskPrimers.primer_field
barcode = params.Mask_Primer_3_MaskPrimers.barcode
revpr = params.Mask_Primer_3_MaskPrimers.revpr
mode = params.Mask_Primer_3_MaskPrimers.mode
failed = params.Mask_Primer_3_MaskPrimers.failed
nproc = params.Mask_Primer_3_MaskPrimers.nproc
maxerror = params.Mask_Primer_3_MaskPrimers.maxerror
umi_length = params.Mask_Primer_3_MaskPrimers.umi_length
start = params.Mask_Primer_3_MaskPrimers.start
extract_length = params.Mask_Primer_3_MaskPrimers.extract_length
maxlen = params.Mask_Primer_3_MaskPrimers.maxlen
skiprc = params.Mask_Primer_3_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_3_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_3_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}else{
	args_1 = args_values[0]
	println args_1
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
	"""
}

}


process collapse_sequences_collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads_unique/$filename"}
input:
 set val(name), file(reads) from g4_11_reads0_g5_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g5_16_reads0_g6_20
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g5_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g5_16_reads_undetermined22
 file "CS_*"  into g5_16_logFile33

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc
failed = params.collapse_sequences_collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process split_sequences_split_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fastq$/) "split_sequence_reads/$filename"}
input:
 set val(name),file(reads) from g5_16_reads0_g6_20

output:
 set val(name), file("*_atleast-*.fastq")  into g6_20_reads0_g7_15

script:
field = params.split_sequences_split_seq.field
num = params.split_sequences_split_seq.num


readArray = reads.toString()

if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num}
"""

}


process Parse_header_table_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "parse_header_table/$filename"}
input:
 set val(name), file(reads) from g6_20_reads0_g7_15

output:
 set val(name),file("*${out}")  into g7_15_reads00

script:
method = params.Parse_header_table_parse_headers.method
act = params.Parse_header_table_parse_headers.act
args = params.Parse_header_table_parse_headers.args


if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	out="_reheader.fastq"
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} ${args}
			"""	
	}else{
		out="_reheader.fastq"
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}


}


process Mask_Primer_3_parse_log_MP {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "MP3_log_tab/$filename"}
input:
 val mate from g_12_mate_g4_9
 set val(name), file(log_file) from g4_11_logFile2_g4_9

output:
 set val(name), file("*.tab")  into g4_9_logFile0_g4_12

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_3_try_report_maskprimer {

input:
 set val(name), file(primers) from g4_9_logFile0_g4_12
 val matee from g_12_mate_g4_12

output:
 file "*.rmd"  into g4_12_rMarkdown0_g4_16


shell:

if(matee=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray.grep(~/.*R1.*/)[0]
	primers_2 = readArray.grep(~/.*R2.*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray.grep(~/.*.tab*/)[0]

	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_3_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "MP3_report/$filename"}
input:
 file rmk from g4_12_rMarkdown0_g4_16

output:
 file "*.html"  into g4_16_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process make_report_pipeline_cat_all_file {

input:
 set val(name), file(log_file) from g0_0_logFile3_g25_0
 set val(name), file(log_file) from g1_0_logFile3_g25_0
 set val(name), file(log_file) from g2_11_logFile3_g25_0
 set val(name), file(log_file) from g3_11_logFile3_g25_0
 set val(name), file(log_file) from g4_11_logFile3_g25_0

output:
 set val(name), file("all_out_file.log")  into g25_0_logFile0_g25_2

script:
readArray = log_file.toString()

"""

echo $readArray
cat out* >> all_out_file.log
"""

}


process make_report_pipeline_report_pipeline {

input:
 set val(name), file(log_files) from g25_0_logFile0_g25_2

output:
 file "*.rmd"  into g25_2_rMarkdown0_g25_1


shell:

readArray = log_files.toString().split(' ')
R1 = readArray[0]

'''
#!/usr/bin/env perl


my $script = <<'EOF';


```{r, message=FALSE, echo=FALSE, results="hide"}
# Setup
library(prestor)
library(knitr)
library(captioner)

plot_titles <- c("Read 1", "Read 2")
if (!exists("tables")) { tables <- captioner(prefix="Table") }
if (!exists("figures")) { figures <- captioner(prefix="Figure") }
tables("count", 
       "The count of reads that passed and failed each processing step.")
figures("steps", 
        paste("The number of reads or read sets retained at each processing step. 
               Shown as raw counts (top) and percentages of input from the previous 
               step (bottom). Steps having more than one column display individual values for", 
              plot_titles[1], "(first column) and", plot_titles[2], "(second column)."))
```

```{r, echo=FALSE}
console_log <- loadConsoleLog(file.path(".","!{R1}"))
```

# Summary of Processing Steps

```{r, echo=FALSE}
count_df <- plotConsoleLog(console_log, sizing="figure")
```

`r figures("steps")`

```{r, echo=FALSE}
kable(count_df[c("step", "task", "total", "pass", "fail")],
      col.names=c("Step", "Task", "Input", "Passed", "Failed"),
      digits=3)
```

`r tables("count")`


EOF
	
open OUT, ">!{name}.rmd";
print OUT $script;
close OUT;

'''
}


process make_report_pipeline_render_rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "pipline_report/$filename"}
input:
 file rmk from g25_2_rMarkdown0_g25_1

output:
 file "*.html"  into g25_1_outputFileHTML00

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
