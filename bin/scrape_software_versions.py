#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'NascentFlow': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'fasterq-dump': ['v_fasterq-dump.txt', r"fasterq-dump : (\S+)"],
    'RseQC': ['v_rseqc.txt', r"RSeQC version(\S+)"],
    'Bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'IGV Tools': ['v_igv-tools.txt', r"IGV Version (\S+)"],
    'BBduk': ['v_bbduk.txt', r"BBDuk Trimming version(\S+)"],
    'Hisat2': ['v_hisat2.txt', r"hisat2-align-s version (\S+)"],
    'preseq': ['v_preseq.txt', r"Preseq version(\S+)"],
    'FStitch': ['v_fstitch.txt', r"FStitch version(\S+)"],
    'Tfit': ['v_tfit.txt', r"Tfit version(\S+)"],
    'DAStk': ['v_dastk.txt', r"DAStk version(\S+)"],
}
    #'RseQC': ['v_rseqc.txt', r"RSeQC version(\S+)"],
    
results = OrderedDict()
results['NascentFlow'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nascentflow-software-versions'
section_name: 'NascentFlow Software Versions'
section_href: 'https://github.com/Dowell-Lab/Nascent-Flow'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
